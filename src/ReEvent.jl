module ReEvent
# Load necessary modules
using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles
using FFTW, DSP, Statistics, BSplineKit
import ColorSchemes.darkrainbow

gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))
cb = darkrainbow

# Include necessary scripts for functions
include("WAVEGEM.jl")
import .WAVEGEM
include("func/text_process.jl"); include("func/signal_processing.jl")
include("func/peak_detect.jl");  include("func/runge_kutta.jl")
include("func/fugees.jl");  include("func/directories.jl")
include("func/wave_theory.jl");  include("func/jonswap.jl")

#############################################################################################
# Global variables
ρ, g, _, d, _, _, γ, fcut = WAVEGEM.GlobInp0
_, case_id, _, _, _, _, _, _, Decpath,_,_,_,_, OFASTpath = WAVEGEM.GlobPaths
T₀ₛ, T₀ₕ, T₀ₚ = WAVEGEM.FOWT[:]
CET, evID = WAVEGEM.CET, WAVEGEM.evID
prop = WAVEGEM.Mprop
frec = WAVEGEM.ReEvflags
Wave = WAVEGEM.Wave

#############################################################################################
# IO directories
if CET == 1        # Fairlead tension critical event
    case_str = "MaxFair"
elseif CET == 2    # Pitch critical event
    case_str = "MaxPitch"
elseif CET == 3    # CoM extreme displacement event 
    case_str = "MaxCoM"
else                # Response to extreme wave event
    case_str = "MaxWave"
end

evdir = joinpath(case_str,"EV$evID")
Eventpath = joinpath(Decpath,evdir)
GRrespath = joinpath(Eventpath,"ERGEEs")
ReEvpath = joinpath(Decpath,"ReEvent")  
fid = joinpath(Eventpath,"event_lin") # Event filename of 1st order elevation
fid_pEWG = joinpath(GRrespath,"EWG_norm_pars")
fid_pG = joinpath(GRrespath,"G_pars")

if !isdir(ReEvpath)
    mkdir(ReEvpath)
end
#############################################################################################
# Read event 1st order elevation
cont = parse_fxw(fid, 0)
tᵢ = (cont[:,1] .- cont[1,1])   # Event time vector
ηᵢ = cont[:,2]      # Surface elevation

## Shift zero time at highest peak
ηHP = findmax(ηᵢ)[1]   # η = η ./ ηHP
iHP = findmax(ηᵢ)[2]
t₀ = tᵢ[iHP]
tᵢ = tᵢ .- t₀

Tₑ = nextpow(2,tᵢ[end]+2*t₀)    # Maximum period / Return period
tₑ = Tₑ    # [s] Simulation duration
Nₛ = Int(tₑ/2)  # No of frequency components
df = 1/tₑ       # Frequency resolution
# fₛ = (Nₛ-1)/tₑ  # [Hz] Sampling frequency
fₛ = 2^3  # [Hz] Sampling frequency

# Read parameters of Gaussian Regression of wave event
cont = parse_fxw(fid_pEWG, 1)   # Read EWGs parameters
t̅ᶜ, T̅ₒ, A̅ₒ, β̃ = cont[:,1], cont[:,2], cont[:,3], cont[:,4]
lenT = length(t̅ᶜ)

cont = parse_fxw(fid_pG, 1)     # Read event and GR parameters
Hₛ, Tₚ, Ω = cont[1:3]

# Read event parameters from file
fpar = joinpath(Eventpath,"ev_pars")
cont = parse_fxw(fpar, 1)  
T₂₋, ΔT = cont[4], cont[5]

# ω̃ = Ω*ωₚⱼ
# T̅ₒ = Tₒ ./ Tₚⱼ;     A̅ₒ = Aₒ ./ Hₛⱼ;     t̅ᶜ = tᶜₒ ./ Tₚⱼ

#############################################################################################
# Propagation of individual EWGs
## Wave characteristic values at peak period Tp
peak_wave = wave_qnts(Tₚ,d)
fₚ = peak_wave.f;   ωₚ = peak_wave.ω    
λₚ = peak_wave.λ;   κₚ = peak_wave.κ
υₚ = peak_wave.υᶜ

## Sampling frequency fs>=2*fcut=fnyq
ωcut = 2π*fcut #  [rad/s] Cut-off frequency for JONSWAP spectrum
Tᵢ = 1/fcut    # [s]
dt = 1/fₛ
Nₜ = Int64(round(tₑ/dt))

## Initialize vectors
t = zeros(Float64,Nₜ)
[t[i] = (i-1)*dt for i ∈ 1:Nₜ]
# t₀ = t[end]/2
# t₀ = 2^6
t = t .- (Tₑ-tᵢ[end]-t₀)

# Re-sampling event
itp_ηᵢ = interpolate(tᵢ, ηᵢ, BSplineOrder(4))
η = itp_ηᵢ.(t)

## Generate surface elevation
TOT = zeros(ComplexF64, Nₜ)  # Total surface elevation
Nₛ₂ = Int(nextpow(2,Nₜ)/2)
SP_TOT = zeros(Float64, Nₛ₂) # Total spectrum of propagated Gaussian EWGs
ωw, ωg = zeros(Float64, lenT), zeros(Float64, lenT)
ωᵍ = 2π/T₂₋
# ωᵍ = π/(t̅ᶜ[1]*Tₚ-t̅ᶜ[2]*Tₚ)
ω₁ = Ω*ωₚ + ωᵍ
ω₂ = Ω*ωₚ - ωᵍ

plt_gn = plot(palette=:darkrainbow)
pltEWG = plot(palette=:darkrainbow)
pltSPEC = plot(palette=:darkrainbow)
plt_demo = plot(palette=[cb[8];cb[11];cb[4]])

for n ∈ 1:lenT
    # Sum and difference frequencies
    ωw[n] = (Ω*ωₚ + 2π/(T̅ₒ[n]*Tₚ))/2
    ωg[n] = (Ω*ωₚ - 2π/(T̅ₒ[n]*Tₚ))/2

    # Propagated EWG
    gₙ, ηᶜ, ηₙ, FR, MAG = recon(prop, A̅ₒ[n], t̅ᶜ[n], T̅ₒ[n], β̃[n], Hₛ, Tₚ, Ω, t)

    TOT[:] = TOT[:] .+ real.(ηₙ)
    SP_TOT[:] = SP_TOT[:] .+ MAG

    plot!(plt_gn, t, gₙ, xlab = "t [s]", ylab = "A [m]", lab = "g$(n)", lw=2, legend=:outerbottom, legendcolumns=min(lenT,8))
    # plot!(plt_gn, xlim=(tᵢ[1]-tₑ/50,tᵢ[end]+tₑ/50))
    # plot!(xlim=(-100,100))
    plot!(pltEWG, t, real(ηₙ), xlab = "t [s]", ylab = "η [m]", lab = "η$(n)", lw=2, legend=:outerbottom, legendcolumns=min(lenT,8))
    # plot!(pltEWG, xlim=(tᵢ[1]-tₑ/50,tᵢ[end]+tₑ/50))
    plot!(pltSPEC, FR, MAG, xlab = L"f [Hz]", ylab = L"A [m]", lab = "g$(n)", lw=2, legend=:outerbottom, legendcolumns=min(lenT,8))
    plot!(pltSPEC, xlim=(0,fcut))

    if n == 1
        plot!(plt_demo, t, gₙ, xlab = "t [s]", ylab = "[m]", lab = L"g_1", lw=2)
        plot!(plt_demo, t, real(ηᶜ), lab="Temporal term", lw=2)
        plot!(plt_demo, t, real(ηₙ), lab="Propagated elementary envelope", lw=2)
        # plot!(plt_demo, xlim=(-4*T̅ₒ[n]*Tₚ,4*T̅ₒ[n]*Tₚ))
    end                         

    fid_EWG = joinpath(GRrespath,"EWG_$n") # File name
    open(fid_EWG, "w")
    # writedlm(fid_EWG, [round.(t*1e6)./1e6 round.(gₙ*1e6)./1e6 round.(real(ηₙ*1e6))./1e6 round.(angle.(ηₙ)*1e6)./1e6], '\t')
    writedlm(fid_EWG, [round.(t*1e6)./1e6 round.(gₙ*1e6)./1e6 round.(ηₙ*1e6)./1e6], '\t')
end
ηᴳ = real(TOT)

# for i ∈ 1:Nₜ
#     if t[i] < -abs(Tᵍ)/4 || t[i] > abs(Tᵍ)/4
#         ηᴳ[i] = 0
#     end
# end

# fr_η, mag_η, ϕ_η, df_η, _ = one_side_asp(ηᵢ, tᵢ)
freq, mag_η, ϕ_η, df_η, _ = one_side_asp(η, t)
_, mag, ϕ, df, _ = one_side_asp(ηᴳ, t)

# Focused Wave (SFWG)
Sfcsd = 1/2 / df .* SP_TOT.^2
A₀ = A̅ₒ[1]*Hₛ

ηfcsd = fcsd_wg(freq, Sfcsd, t, A₀, 0, 0)
_, mag_fcsd, pfi, _ = one_side_asp(ηfcsd, t)

# NewWave
fⱼ,Sⱼ = spectrum(Hₛ,Tₚ,γ,1/fcut,Tₑ,Nₛ) # Generate JONSWAP spectrum
ηⱼ = fcsd_wg(fⱼ, Sⱼ, t, Hₛ, 0, 0)

# BEAT
ηᶜ,Tᵍ,Tʷ = wave_group(1,2π/ω₁,1,2π/ω₂,d,t,0, 0)
# ηᶜ,_ = wave_group(1,Tₚ,1,2π/(Ω*ωₚ),d,t,0, 0)
ηᴮ = A̅ₒ[1]*Hₛ .* ηᶜ

for i ∈ 1:Nₜ
    if t[i] < -1/4*abs(Tᵍ) || t[i] > 1/4*abs(Tᵍ)
        ηᴮ[i] = 0
    end
end
_, magᴮ, ϕᴮ, _ = one_side_asp(ηᴮ, t)

# DWG
p⁺ = findmax(A̅ₒ.*T̅ₒ)[2]
gₙ = elem_wg(t, A̅ₒ*Hₛ, zeros(Float64,lenT), T̅ₒ*Tₚ, p⁺)
# G = gauss_fun(t, A̅ₒ*Hₛ, t̅ᶜ*Tₚ, T̅ₒ*Tₚ);
# ηᶜ = exp.(-1im * Ω*ωₚ * t)
ηᶜ,Tᵍ,Tʷ = wave_group(1,2π/ω₁,1,2π/ω₂,d,t,t̅ᶜ[p⁺]*Tₚ, 0)
ηₙ = gₙ.*real(ηᶜ)
# ηₙ = G.*real(ηᶜ)
_, magD, ang, df,_ = one_side_asp(ηₙ,t)
Sᵥ = 0.5/df .* magD.^2
ηᴰ = fcsd_wg(freq, Sᵥ, t, A̅ₒ[p⁺]*Hₛ, 0, 0)

# Statistical measures
## Mean Square Errors
MSEₜ = sum((η.-ηᴳ).^2) / Nₜ 
MSEₛ = sum((mag_η.-mag).^2) / Nₛ₂
## Cross-correlation
CCRₜ = maximum(xcorr(η, ηᴳ))
CCRₛ = maximum(xcorr(mag_η, mag))

## Cosine Similarity
Na = sqrt(sum(η.^2));   Nb = sqrt(sum(ηᴳ.^2))
CSₜ = dot(η, ηᴳ) / (Na*Nb) 
Na = sqrt(sum(mag_η.^2));   Nb = sqrt(sum(mag.^2))
CSₛ = dot(mag_η, mag) / (Na*Nb) 

println("MSEₜ = $MSEₜ")
println("CCRₜ = $CCRₜ")
println("CSₜ = $CSₜ")
println("MSEₛ = $MSEₛ")
println("CCRₛ = $CCRₛ")
println("CSₛ = $CSₛ")

# Plot original and reconstructed signals
plt_eta = plot(xlab=L"t~[s]", ylab=L"\eta ~[m]", palette=[cb[8];cb[11];cb[4]]) 
# plot!(tᵢ,ηᵢ, lab = L"\eta(t)",lw=2)
plot!(t,η, lab = L"\eta(t)",lw=2)
plot!(t, ηᴳ, lab = L"\eta^{\prime}_G (t)", lw=2, line=:dot)
# plot!(xlim=(tᵢ[1]-tₑ/50,tᵢ[end]+tₑ/50))
# plot!(xlim=(-100,100))

plt_fcsd = plot(xlab=L"t~[s]", ylab=L"\eta ~[m]", palette=[cb[8];cb[11];cb[4];cb[1];cb[6]]) 
plot!(t,η, lab = L"\eta(t)",lw=2)
plot!(t, ηfcsd, lab = L"\eta_{fcsd}", lw=2, line=:dot)
plot!(t, ηᴰ, lab = L"\eta_D", lw=2, line=:dot)
plot!(t, ηᴮ, lab = L"\eta_B", lw=1, line=:dash)
plot!(t, ηⱼ, lab = L"\eta_{NW}", lw=1)
# plot!(xlim=(tᵢ[1]-tₑ/50,tᵢ[end]+tₑ/50))
# plot!(xlim=(-100,100))

# Spectra of resulting signals
plt_spec = plot(xlab = L"f~[Hz]", ylab = L"S~ [m]", palette=[cb[8];cb[11];cb[4];cb[1]], legend=:topright)
plot!(freq, mag_η,  lab = L"S_{ev}", lw=2)
plot!(freq, mag, lab = L"S_G(f)", lw=2, line=:dot)
plot!(freq, SP_TOT, lab = L"S(f)", lw=2, line=:dashdot)
# plot!(freq, magD, lab = L"S_D", lw=2, line=:dash)
plot!(freq, magᴮ, lab = L"S_B", lw=2, line=:dash)
# plot!(freq, mag_fcsd, lab = "Focused", lw=2, line=:dashdot)
plot!(xlim=(0,0.625))

# ifft
H2 = Nₛ₂*mag.*exp.(1im*ϕ)
H2 = vcat(H2,0,conj(H2[end:-1:2]))
ETA = real(ifft(H2))

plt_ifft = plot(xlab=L"t~[s]", ylab=L"\eta ~[m]", palette=[cb[11];cb[8]]) 
plot!(t, real(TOT), lab = L"\sum EWGs", lw=2, line=:dot)
plot!(t, ETA, lab="ifft")
# plot!(xlim=(tᵢ[1]-tₑ/100,tᵢ[end]+tₑ/100))

display(plt_gn)
display(pltEWG)
display(pltSPEC)
display(plt_eta)
display(plt_fcsd)
display(plt_spec)
display(plt_demo)
display(plt_ifft)

if frec 
    make_dirs(2, joinpath(ReEvpath,case_str), joinpath(ReEvpath,evdir))

    fid = joinpath(ReEvpath,evdir,"reconstruction")
    savefig(plt_eta, fid*".svg");  savefig(plt_eta, fid*".png")
    fid = joinpath(ReEvpath,evdir,"focused_wave")
    savefig(plt_fcsd, fid*".svg");  savefig(plt_fcsd, fid*".png")
    fid = joinpath(ReEvpath,evdir,"spectra")
    savefig(plt_spec, fid*".svg");  savefig(plt_spec, fid*".png")
    
    fid = joinpath(ReEvpath,evdir,"elem_envs")
    savefig(plt_gn, fid*".svg");  savefig(plt_gn, fid*".png")
    
    fid = joinpath(ReEvpath,evdir,"EWGs")
    savefig(pltEWG, fid*".svg");  savefig(pltEWG, fid*".png")

    fid = joinpath(ReEvpath,evdir,"event.Elev")
    open(fid, "w")
    writedlm(fid, [t.+t₀ η], '\t')

    suffix = "resamp_FR.fronts"
    fid = joinpath(ReEvpath,evdir,suffix)
    open(fid, "w")
    head0 = [suffix "" "" ""]
    head1 = ["f" "a" "angle" "phase"]
    head2 = ["Hz" "m" "rad" "rad"]
    wcont = [head0; head1; head2; freq mag_η zeros(Float64,Nₛ₂) ϕ_η]
    writedlm(fid, wcont, '\t')

    # OpenFAST elevation input files
    OFElev = joinpath(OFASTpath,"ExtElev","$(case_id)","$(WAVEGEM.run_id)")
    make_dirs(2, joinpath(OFElev,case_str), joinpath(OFElev,evdir))
    
    fid = joinpath(OFElev,evdir,"event.Elev")
    open(fid, "w")
    writedlm(fid, [t.+(Tₑ-tᵢ[end]-t₀) η], '\t')

    fid = joinpath(OFElev,evdir,"ReFoGWs.Elev")
    open(fid, "w")
    writedlm(fid, [t.+(Tₑ-tᵢ[end]-t₀) ηᴳ], '\t')

    suffix = "ReFoGWs_FR.fronts"
    fid = joinpath(ReEvpath,evdir,suffix)
    open(fid, "w")
    head0 = [suffix "" "" ""]
    head1 = ["f" "a" "angle" "phase"]
    head2 = ["Hz" "m" "rad" "rad"]
    wcont = [head0; head1; head2; freq mag zeros(Float64,Nₛ₂) ϕ]
    writedlm(fid, wcont, '\t')

    fid = joinpath(OFElev,evdir,"SFWG.Elev")
    open(fid, "w")
    writedlm(fid, [t.+(Tₑ-tᵢ[end]-t₀) ηfcsd], '\t')

    suffix = "FCSD_FR.fronts"
    fid = joinpath(ReEvpath,evdir,suffix)
    open(fid, "w")
    head0 = [suffix "" "" ""]
    head1 = ["f" "a" "angle" "phase"]
    head2 = ["Hz" "m" "rad" "rad"]
    wcont = [head0; head1; head2; freq mag_fcsd zeros(Float64,Nₛ₂) zeros(Float64,Nₛ₂)]
    writedlm(fid, wcont, '\t')

    fid = joinpath(OFElev,evdir,"NW.Elev")
    open(fid, "w")
    writedlm(fid, [t.+(Tₑ-tᵢ[end]-t₀) ηⱼ], '\t')

    fid = joinpath(OFElev,evdir,"BEAT.Elev")
    open(fid, "w")
    writedlm(fid, [t.+(Tₑ-tᵢ[end]-t₀) ηᴮ], '\t')

    fid = joinpath(OFElev,evdir,"DWG.Elev")
    open(fid, "w")
    writedlm(fid, [t.+(Tₑ-tᵢ[end]-t₀) ηᴰ], '\t')
end

end
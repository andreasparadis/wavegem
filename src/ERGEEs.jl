module ERGEEs 
# Load necessary modules
using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles
using FFTW.AbstractFFTs,  DSP, BSplineKit
using Statistics: std, mean
using CurveFit: linear_fit

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
ρ, g, _, d, _, _, _, fcut = WAVEGEM.GlobInp0
T₀ₛ, T₀ₕ, T₀ₚ = WAVEGEM.FOWT[:]

#############################################################################################
# IO directories
Parent, _, _, _, _, _, _, _, Decpath, _ = WAVEGEM.GlobPaths

Eventpath = joinpath(Decpath,WAVEGEM.evdir)
fid_res = joinpath(Eventpath,"ERGEEs")  
fid_1st = joinpath(Eventpath,"event_lin") # Path to event file of 1st order elevation
fpar = joinpath(Decpath,"sim_pars")       # Simulation parameters file path

make_dirs(1, fid_res)

#############################################################################################
# Module specific inputs
## Read simulation parameters from file
cont = parse_fxw(fpar, 1)  
Hₛ, Tₚ = cont[1], cont[2]
Tcut = 1/fcut

frec, ftrunc = WAVEGEM.ERGEEflags 
fcˡᵖ, std_fac, MinPeakDist, ϵ, Nτ, dτ, Tlb, Tub = WAVEGEM.ERGEEsInp

#############################################################################################
# Pre-processing of wave event
## Read event surface elevation from file
cont = parse_fxw(fid_1st, 0)

tᵢ = (cont[:,1] .- cont[1,1]) # initial time vector
ηᵢ = cont[:,2] # Surface elevation
Lᵢ = length(tᵢ)

## Shift zero time at highest peak
ηHP = findmax(ηᵢ)[1]   # η = η ./ ηHP
iHP = findmax(ηᵢ)[2]
t = tᵢ .- tᵢ[iHP]     

# Truncate event based on number of up-crossings before and after highest peak
iDᵤ = findall(diff(sign.(ηᵢ)) .== 2) # Up-crossings of surface elevation
LiDᵤ = length(iDᵤ)

istart, iend = Int(0), Int(0)
cnt = Int(iHP - minimum(abs.(iHP.-iDᵤ)))

if ftrunc
    # Select based on No of up-crossings before and after
    for i ∈ 1:LiDᵤ
        if iDᵤ[i] == cnt
            icnt = i
            global istart = iDᵤ[icnt-WAVEGEM.pNuc]
            global iend = iDᵤ[icnt+WAVEGEM.aNuc]
        end
    end
else
    # From 1st to last upcrossing
    global istart = iDᵤ[1]
    global iend = iDᵤ[LiDᵤ]
end

## Truncate event
t = t[istart:iend]
η = ηᵢ[istart:iend]
M = length(t);      Øᴹ = zeros(Float64,M) 
dt = (t[end]-t[1])/(M-1)

## Event spectral analysis and peak frequency
fr_η, mag_η, ϕ_η,_ = one_side_asp(η, t)
fᵉₚ = fr_η[findmax(mag_η)[2]]

#############################################################################################
# Hilbert transform and envelope of surface elevation
𝓗 = hilbert(η);        uᵢ = abs.(𝓗)

# Low-pass filtering of envelope
fₛˡᵖ = round(1/dt)
Ntaps = nextpow(2,fₛˡᵖ/fcˡᵖ)   # No of taps
u = low_pass_filter(uᵢ,fcˡᵖ,fₛˡᵖ,Ntaps)

# Peak detection based on 2nd derivative of signal
MinPeakVal = std_fac*std(η) 
A, tᶜ, i⁺, uₜ, uₜₜ = peaks_max_ext(u,t, MinPeakVal, MinPeakDist,true)

#############################################################################################
# Up- and down-crossings
iDue = findall(diff(sign.(u.-MinPeakVal)) .== 2) # Up-crossings of envelope
iDde = findall(diff(sign.(u.-MinPeakVal)) .== -2) # Down-crossings of envelope

TIE, TDE, Tz = [0], [0], 0
if length(iDue) > 1
    TIE = [abs(t[iDue[i]] - t[iDue[i+1]]) for i ∈ 1:length(iDue)-1]
    TDE = [abs(t[iDue[i]] - t[iDde[i]]) for i ∈ 1:length(iDue)]
    Tz = mean(TIE)
else
    Tz = TIE
end
TD = (t[end]-t[1])
#############################################################################################
# Perform the least squares fitting
## Initial condition for L based on ∂ₜ²g(tᶜ) = ∂ₜ²u(tᶜ)
N = length(tᶜ)     
Øᴺ = zeros(Float64,N);  Ø = Array{Float64}(undef,0)
T = Øᴺ

for n ∈ 1:N
    i = i⁺[n]
    if uₜₜ[i] < 0 
        T[n] = sqrt.(-2*A[n]./uₜₜ[i])
    else
        T[n] = sqrt.(2*A[n]./uₜₜ[i])
    end
end

# Definition of the ODE's dLₘ 
τ = zeros(Float64,Nτ)
[τ[i] = (i-1)*dτ for i ∈ 2:Nτ]
println("τₑ = ", τ[end])
fun(τ,T) = dLdτ(τ,T, t, A, tᶜ,u, dt)

# Solve optimization problem
Tₒₚₜ = RK4_nln_sys(τ,T,fun,Tlb,Tub,ϵ)

# Optimal values
Tₒ, Aₒ, tᶜₒ = [Ø[:] for _ = 1:3]
i⁺ₒ = Array{Int64}(undef,0)
pltCONV = plot(palette=:darkrainbow)

for n ∈ 1:N
    plot!(pltCONV, Tₒₚₜ[2][:,n], xlab = "Iteration",  ylab = L"T_o^i [s]", xscale=:log10, lab = "To[$(n)]", lw=2, legend=:outerbottom, legendcolumns=N)
    if Tₒₚₜ[1][n] > Tcut
        push!(Tₒ, Tₒₚₜ[1][n])
        push!(Aₒ, A[n])
        push!(tᶜₒ, tᶜ[n])
        push!(i⁺ₒ,i⁺[n])
    end
end

lenT = length(Tₒ)

#############################################################################################
# Resulting Gaussian approximation of the envelope
G = gauss_fun(t, Aₒ, tᶜₒ, Tₒ)

# Propagation of resulting WGs
θ = angle.(𝓗)       # θ = ϕ-ωt
eθ = exp.(1im*θ)     # exp(iθ) = exp(-iωt)*exp(iϕ)
# ηₜ = G.*eθ

## Instantaneous frequency
omH = zeros(Float64, M)
θᵤ = unwrap(θ)
for i ∈ 2:M-1
    omH[i] = (θᵤ[i+1] - θᵤ[i-1]) / (t[i+1]-t[i-1])
end

frH, mgH, angH, _,_ = one_side_asp(omH,t)

ωₒ = 2π ./ Tₒ
ωᵢ = zeros(Float64, lenT)
βᵢ = zeros(Float64, lenT)
βₒ = zeros(Float64, lenT)

for n ∈ 1:lenT
    eϕᵢ = eθ .* exp.(1im*ωₒ[n]*t)           # exp(iϕ) = exp(iθ)/exp(-iωt) = exp(iθ)*exp(iωt)
    ϕᵢ = real(-1im*log.(eϕᵢ))
    ϕᵢᵢ = unwrap(ϕᵢ)
    # ωᵢ[n] = (ϕᵢᵢ[end]-ϕᵢᵢ[1])/(t[end]-t[1])
    local fit = linear_fit(t,ϕᵢᵢ)
    ωᵢ[n] = fit[2]    # ϕ = ω⁺*t + β
    βᵢ[n] = fit[1]

    nᵦ = Int(round((tᶜₒ[n]-t[1])/dt + 1))
    βₒ[n] = ϕᵢᵢ[nᵦ]
    # βₒ[n] = ωᵢ[n]*tᶜₒ[n]+βᵢ[n]
end

ω̃ = (ωᵢ.-ωₒ)[1]
β̃ = (βᵢ.-βₒ)

eωt = sum(exp.(-1im*ωₒ[n]*t) for n ∈ 1:lenT)    # exp(-iωt) = ∑exp(-iωₙt)
# eωt = eωt ./ maximum(real(eωt))
ωt = unwrap(real(-1im*log.(eωt)))
f⁺ᵢ = ωᵢ/2/π
T⁺ᵢ = 1/f⁺ᵢ

eϕ = eθ ./ eωt
ϕ = real(-1im*log.(eϕ))
fit = linear_fit(t,unwrap(ϕ))
ω⁺ = fit[2]    # ϕ = ω⁺*t + β ≡ mean(ωᵢ)
β = fit[1]     # ≡ mean(βᵢ)
f⁺ = ω⁺/2/π
T⁺ = 1/f⁺

eϕ̃ = exp.(1im*(ω⁺*t .+ β))
eθ̃ = exp.(1im*(ω̃*t .+ β̃[1]))
# eθ̃ = exp.(1im*(ω⁺*t .+ β)) .* exp.(1im*ωt)
# eθ̃ = exp.(1im*(ω⁺*t .+ β)) .* eωt
ηₜ = G.*real(eθ̃)

#############################################################################################
# Normalization
ωₚ = 2π/Tₚ
Ω = ω̃/ωₚ;      ω̅⁺ = ω⁺/ωₚ;      ω̅ᵢ = ωᵢ./ωₚ
G̅ = G ./ maximum(u)
t̅ = t ./ Tₚ
T̅z = Tz / Tₚ    
T̅IE = TIE ./ Tₚ;   T̅DE = TDE ./ Tₚ
T̅ₒ = Tₒ ./ Tₚ;     A̅ₒ = Aₒ ./ Hₛ;     t̅ᶜ = tᶜₒ ./ Tₚ

t̅ₘ = zeros(Float64, lenT)   # Peak times based on envelope intervals (from up-crossings)
t̅ₘ[1] = -T̅IE[1]
if length(TIE) > 1
    for n ∈ 1:length(TIE)-1
        t̅ₘ[n+1] = t̅ₘ[n] + T̅IE[n]
    end
end

#############################################################################################
# Propagation of individual EWGs
TOT = zeros(ComplexF64, M)  # Total surface elevation
Nₛ = Int(nextpow(2,M)/2)
SP_TOT = zeros(Float64, Nₛ) # Total spectrum of propagated Gaussian EWGs

plt_gn = plot(palette=:darkrainbow)
pltEWG = plot(palette=:darkrainbow)
pltSPEC = plot(palette=:darkrainbow)

for n ∈ 1:lenT  
    # Propagated EWG - Focused WG
    gₙ, ηᶜ, ηₙ, FR, MAG = recon(0, A̅ₒ[n], t̅ᶜ[n], T̅ₒ[n], β̃[n], Hₛ, Tₚ, Ω, t)

    TOT[:] = TOT[:] .+ ηₙ
    SP_TOT[:] = SP_TOT[:] .+ MAG

    plot!(plt_gn, t̅*Tₚ, gₙ, xlab = "t [s]", ylab = "A [m]", lab = "g$(n)", lw=2, legend=:outerbottom, legendcolumns=lenT)
    plot!(pltEWG, t̅*Tₚ, real(ηₙ), xlab = "t [s]", ylab = "η [m]", lab = "η$(n)", lw=2, legend=:outerbottom, legendcolumns=lenT)
    plot!(pltSPEC, FR, MAG, xlab = L"f [Hz]", ylab = L"A [m^2 s]", lab = "g$(n)", lw=2, legend=:outerbottom, legendcolumns=lenT)                            

    fid_EWG = joinpath(fid_res,"EWG_$n") # File name
    open(fid_EWG, "w")
    # writedlm(fid_EWG, [round.(t*1e6)./1e6 round.(gₙ*1e6)./1e6 round.(real(ηₙ*1e6))./1e6 round.(angle.(ηₙ)*1e6)./1e6], '\t')
    writedlm(fid_EWG, [round.(t*1e6)./1e6 round.(gₙ*1e6)./1e6 round.(ηₙ*1e6)./1e6], '\t')
end

#############################################################################################
# PLOTS
# Signal and initial envelope
plt_sigenv = plot(xlab=L"t~[s]", ylab=L"[m]", palette=[cb[11];cb[4];cb[8]])
plot!(t, [η uᵢ], lab=[L"η(t)" L"u(t)"], lw=[2 2])
plot!(t, MinPeakVal*ones(Int,M), lab="$std_fac std(η)", line=:dash)

# Filtered vs initial envelope
plt_env = plot(xlab=L"t~[s]", ylab=L"[m]", palette=[cb[11];cb[4];cb[8]])
plot!(t, [uᵢ u], lab=[L"u(t)" L"u_{flt}(t)"], lw=[2 2], line=[:solid :dashdot])
# plot!(t, [uᵢ uᵢᵢ u], lab=[L"u(t)" L"u_{ii}(t)" L"u_{flt}(t)"], lw=[2 2 2], line=[:solid :solid :solid])
# plot!(t, std_fac_u*std(uᵢᵢ)*ones(Int,M), lab="$std_fac_u std(u)", line=:dash)
plot!(t, MinPeakVal*ones(Int,M), lab="$std_fac std(η)", line=:dash)
# plot!(t, MinPeakVal*ones(Int,M), lab=L"H_s/2", line=:dash)

# Gaussian approximation against envelope and scatter plot of peaks
plt_G = plot(xlab = L"t~[s]", ylab = L"[m]", palette=[cb[4];cb[8];cb[11]], legend=:outertop, legendcolumns=4)
plot!(t, [u real(G)], lab = [L"u_{flt}(t)" L"G(t)"], lw=[2 2])
plot!(t, MinPeakVal*ones(Int,M), lab="$std_fac std(η)", line=:dash)
# plot!(t, MinPeakVal*ones(Int,M), lab=L"H_s/2", line=:dash)
plot!(tᶜ, A, seriestype=:scatter, ms=2, mc=:red, lab = "Peaks")

# Plot EWGs
# plot!(pltEWG, t, η, lab = L"\eta(t)",lw=2, line=:dot)
plot!(pltSPEC, xlim=(0,fcut))

# Plot original and reconstructed signals
plt_eta0 = plot(xlab=L"t~[s]", ylab=L"\eta ~[m]", palette=[cb[11];cb[4]]) 
plot!(t, η, lab = L"\eta(t)",lw=2)
plot!(t, real(ηₜ), lab=L"\eta_G(t)", lw=2, line=:dot)

plt_eta1 = plot(xlab=L"t~[s]", ylab=L"\eta ~[m]", palette=[cb[11];cb[8]]) 
plot!(t,η, lab = L"\eta(t)",lw=2)
plot!(t̅*Tₚ, real(TOT), lab = L"\sum EWGs", lw=2, line=:dot)

plt_eta2 = plot(xlab=L"t~[s]", ylab=L"\eta ~[m]", palette=[cb[4];cb[8]]) 
plot!(t, real(ηₜ), lab=L"\eta_G(t)", lw=2, line=:dot)
plot!(t̅*Tₚ, real(TOT), lab = L"\sum EWGs", lw=2)

# Spectra of resulting signals
plt_spec = plot(palette=[cb[11];cb[4];cb[8];cb[1]], legend=:topright)
plot!(fr_η, mag_η, xlab = "f [Hz]", ylab = L"A [m^2 s]", lab = "Original", lw=2)
freq, mag, _ = one_side_asp(real(ηₜ), t)
plot!(freq, mag, xlab = "f [Hz]", ylab = L"A [m^2 s]", lab = L"\eta_{new}", lw=2)
freq, mag, _ = one_side_asp(real(TOT), t̅*Tₚ)
plot!(freq, mag, xlab = "f [Hz]", ylab = L"A [m^2 s]", lab = "Sum", lw=2, line=:dot)
plot!(freq, SP_TOT, xlab = "f [Hz]", ylab = L"A [m^2 s]", lab = "SP_TOT", lw=2, line=:dashdot)
plot!(xlim=(0,fcut))

# Plot the group parameters in ascending order
plt_pars = plot(xlab = L"T_o~[s]", ylab = L"A_o~[m]")
plot!(Tₒ, Aₒ, seriestype=:scatter, ms=2, mc=:red, lab = "Group parameters (sorted)")

# Normalized envelope
plt_G̅ = plot(t̅, G̅)
plot!(plt_G̅, xlab=L"t/T_p~[-]", ylab=L"\overline{G} = G/u_{max}~[-]")

plt_ang_un = plot(xlab=L"t/T_p", ylab=L"\angle ~[rad]", palette=[cb[11];cb[4];cb[8]])
plot!(plt_ang_un, t̅, θᵤ, lab=L"\theta (t)", lw=2)
plot!(plt_ang_un, t̅, ω⁺*t .+ β .+ ωt, lab=L"\theta_{cal} (t)", lw=2)
plot!(plt_ang_un, t̅, ω⁺*t .+ β, lab=L"\phi_{cal} (t)", line=:dash)
plot!(plt_ang_un, t̅, -ωt, lab=L"ωt_{cal} (t)", line=:dash)

plt_ang = plot(xlab=L"t/T_p", ylab="Temporal term ~[-]", palette=[cb[11];cb[4];cb[8]])
plot!(plt_ang, t̅, real(eθ), lab=L"Re\{e^{\theta (t)}\}", lw=2)
plot!(plt_ang, t̅, real(eθ̃), lab=L"Re\{e^{\theta_{cal} (t)}\}", lw=2)
plot!(plt_ang, t̅, real(eϕ̃), lab=L"Re\{e^{\phi_{cal} (t)}\}", line=:dash)
plot!(plt_ang, t̅, real(exp.(1im*ωt)), lab=L"Re\{e^{ωt_{cal} (t)}\}", line=:dash)

plot!(pltCONV,xlim=(1e1,1e5))

# Display plots
display(plt_sigenv)
display(plt_env)
display(plt_G)
display(plt_gn)
display(pltEWG)
display(pltSPEC)
display(plt_eta0)
display(plt_eta1)
display(plt_eta2)
display(plt_spec)
display(plt_pars)
display(plt_G̅)
display(plt_ang_un)
display(plt_ang)
display(pltCONV)

etaG = zeros(Float64, Lᵢ)
etaG[istart:iend] = real(ηₜ)[:]
plot(tᵢ, ηᵢ, lw=2, lab=L"\eta (t)")
# plot!(tᵢ,etaG, lw=2, lab=L"G (t)")

EWGsum = zeros(Float64, Lᵢ)
EWGsum[istart:iend] = real(TOT)[:]
display(plot!(tᵢ,EWGsum, lw=2, lab=L"\sum g_i (t)"))

ηTOT = fcsd_wg(freq, SP_TOT, t, Aₒ[1], 0, β̃[1])
etaFCSD = zeros(Float64, Lᵢ)
etaFCSD[istart:iend] = ηTOT[:]
# display(plot!(tᵢ, etaFCSD, lw=2, lab="Focused wave"))

# Save figures
fid = joinpath(fid_res,"envelope")
savefig(plt_sigenv, fid*".svg");  savefig(plt_sigenv, fid*".png")
fid = joinpath(fid_res,"env_Ginterp")
savefig(plt_G, fid*".svg");  savefig(plt_G, fid*".png")
fid = joinpath(fid_res,"EWGs")
savefig(pltEWG, fid*".svg");    savefig(pltEWG, fid*".png")
fid = joinpath(fid_res,"new_sig_comp")
savefig(plt_eta0, fid*".svg");   savefig(plt_eta0, fid*".png")
fid = joinpath(fid_res,"EWGpars")
savefig(plt_pars, fid*".svg");   savefig(plt_pars, fid*".png")

#############################################################################################
if frec
    # Write results in text files
    t = round.(t*1e6)./1e6
    ηₜ = round.(ηₜ*1e6)./1e6
    TOT = round.(TOT*1e6)./1e6
    tᶜₒ = round.(tᶜₒ*1e6)./1e6
    Aₒ = round.(Aₒ*1e6)./1e6
    Tₒ = round.(Tₒ*1e6)./1e6

    fid = joinpath(fid_res,"etaG")
    open(fid, "w")
    writedlm(fid, [tᵢ etaG], '\t')
    # writedlm(fid, [t real(ηₜ) angle.(ηₜ)], '\t')

    fid = joinpath(fid_res,"EWGsum")
    open(fid, "w")
    writedlm(fid, [tᵢ EWGsum], '\t')
    # writedlm(fid, [t real(TOT) angle.(TOT)], '\t')

    fid = joinpath(fid_res,"etaFCSD")
    open(fid, "w")
    writedlm(fid, [tᵢ etaFCSD], '\t')

    fid = joinpath(fid_res,"EWGpars")
    head = ["tᶜₒ [s]" "Tₒ [s]" "Aₒ [m]"]
    open(fid, "w")
    writedlm(fid, [head; tᶜₒ Tₒ Aₒ], '\t')    

    fid = joinpath(fid_res,"EWG_norm_pars")
    # head = ["t̅ᶜ [-]" "T̅ₒ [-]" "A̅ₒ [-]" "βᵢ [rad]" "ω̅ᵢ [-]" "T̅DE [-]" "t̅ₘ"]
    head = ["t̅ᶜ [-]" "T̅ₒ [-]" "A̅ₒ [-]" "β̃ [rad]" "ω̅ᵢ [-]" "T̅ₛ" "T̅ₕ" "T̅ₚ"]
    open(fid, "w")
    # row = round.([t̅ᶜ T̅ₒ A̅ₒ β̃ ω̅ᵢ T̅DE t̅ₘ]*1e6)./1e6
    row = round.([t̅ᶜ T̅ₒ A̅ₒ β̃ ω̅ᵢ Tₒ./T₀ₛ Tₒ./T₀ₕ Tₒ./T₀ₚ]*1e6)./1e6
    writedlm(fid, [head; row], '\t')

    fid = joinpath(fid_res,"G_pars")
    head = ["Hₛ [m]" "Tₚ [s]" "Ω [-]" "β [rad]" "$std_fac*σₙ" "ηₘₐₓ [m]" "uₘₐₓ [m]" "t̅ₑᵥ [-]" "ω̅⁺ [-]" "T̅z [-]"]
    open(fid, "w")
    row = round.([Hₛ Tₚ Ω β MinPeakVal ηHP maximum(u) t̅[end]-t̅[1] ω̅⁺ T̅z]*1e6)./1e6
    writedlm(fid, [head; row], '\t')
end

end
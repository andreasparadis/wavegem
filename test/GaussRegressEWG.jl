# Load necessary modules
using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles
using FFTW, DSP, LsqFit, Statistics,  Random
using Distributions: Rayleigh, Uniform, Normal
using CurveFit: curve_fit, linear_fit, Polynomial, expsum_fit, ExpFit
import ColorSchemes.darkrainbow

gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))
cb = darkrainbow

# Include necessary scripts for functions
include("text_process.jl")
include("signal_processing.jl")
include("peak_detect.jl")
include("runge_kutta.jl")
include("Gauss_Reg_t.jl")
include("directories.jl")
include("wave_theory.jl")

#############################################################################################
# Global variables
Hₛⱼ = 7.306;   Tₚⱼ = 9.9
# Hₛⱼ = 0.069;  Tₚⱼ = 1       
ωₚⱼ = 2π/Tₚⱼ

# Semi-submersible FOWT eigenfrequencies
T₀ₛ¹ = 250.0000;    T₀ₛ² = 35.714107  
T₀ₕ = 17.543772;    T₀ₚ = 27.026892

g, d ::Float64 = 9.81, 100

# IO directories
dir0::String = "library/SE/LCH/1/"
fid0::String = dir0*"Decomposition/MaxWave/"
fid_res::String = fid0*"GR_Results/"
fid_ev::String = "event_lin" # File name

# dir0::String = "/home/andreasp/WAVEGEM/library/JFM/"
# fid0::String = dir0*"EWG/"                      
# fid_res::String = fid0*"GR_Results/"
# fid_ev::String = "2_1.txt" # File name    

# dir0::String = "/home/andreasp/WAVEGEM/library/JFM/EV_2/"
# fid0::String = dir0*"TGTs/"                      
# fid_res::String = fid0*"GR_Results/"
# fid_ev::String = "2_G1.txt" # File name   

if !isdir(fid_res)
    mkdir(fid_res)
end
#############################################################################################
# Pre-processing of wave event
## Read event surface elevation from file
fid = fid0*fid_ev # Path to file
cont = parse_fxw(fid, 0)

tᵢ = round.((cont[:,1] .- cont[1,1])*100)/100 # x-coordinates vector
ηᵢ = cont[:,2] # Surface elevation
Lᵢ = length(tᵢ)

## Shift zero time at highest peak
ηHP = findmax(ηᵢ)[1]   # η = η ./ ηHP
iHP = findmax(ηᵢ)[2]
t = tᵢ .- tᵢ[iHP]     

# Truncate event based on number of up-crossings before and after highest peak
iDᵤ = findall(diff(sign.(ηᵢ)) .== 2) # Up-crossings of surface elevation

## Select No of up-crossings before and after
istart, iend = 0, 0
cnt = iHP - minimum(abs.(iHP.-iDᵤ))

for i ∈ 1:length(iDᵤ) 
    if iDᵤ[i] == cnt
        icnt = i
        global istart = iDᵤ[icnt-14]
        global iend = iDᵤ[icnt+5]
    end
end

## Truncate event
t = t[istart:iend]
η = ηᵢ[istart:iend]
M = length(t);      Øᴹ = zeros(Float64,M) 
dt = (t[end]-t[1])/(M-1)

## Event spectral analysis and peak frequency
fr_η, mag_η, _ = one_side_asp(η, t)
fₚ = fr_η[findmax(mag_η)[2]];   Tₚ = 1/fₚ;  ωₚ = 2π*fₚ

# fcut = round(4*fₚ*10)/10 # Cut-off frequency (=4fₚ) [Hz]
fcut = 0.625
Tcut = 1/fcut

#############################################################################################
# Hilbert transform and envelope of surface elevation
𝓗 = hilbert(η);        uᵢ = abs.(𝓗)
## Define filter cut-off frequency based on peak frequency of envelope
# fr_uᵢ, mag_uᵢ = one_side_asp(uᵢ, t)
# iEP = findmax(mag_uᵢ)[2]
# fᴴcut = 2*fr_uᵢ[iEP]
fᴴcut = fₚ

# Low-pass filtering of envelope
fₛˡᵖ = 64*fcut; Ntaps = 2^9+1   # Sampling frequency and No of taps of low-pass filter
u = low_pass_filter(uᵢ,fcut,fₛˡᵖ,Ntaps)

# Peak detection based on 2nd derivative of signal - Define A
std_fac = 1.4 # Coefficient of standard deviation
MinPeakVal, MinPeakDist = std_fac*std(η), 7*Tcut
# MinPeakVal, MinPeakDist = 2*Hₛⱼ/4, 2*Tcut
A, tᶜ, i⁺, uₜ, uₜₜ = peaks_max_ext(u,t, MinPeakVal, MinPeakDist)

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
N = length(tᶜ);     Øᴺ = zeros(Float64,N)
T = Øᴺ

for n ∈ 1:N
    i = i⁺[n]
    if uₜₜ[i] < 0 
        T[n] = sqrt.(-2*A[n]./uₜₜ[i])
    else
        T[n] = sqrt.(2*A[n]./uₜₜ[i])
    end
end
# T = 2π/ωₚⱼ*ones(Float64,N)
# tᶜ = tᶜ .+ [-0.5;-1]
# tᶜ[4] = (tᶜ[4] + tᶜ[7])/2
# T[4] = 0.0001
# A[4] = (A[4]+A[7])/2
# A[4] = 0

# Definition of the ODE's dLₘ 
ϵ = 1e-5
Nτ = Int(1e5)
dτ = 1000*ϵ 
τ = zeros(Float64,Nτ)
[τ[i] = (i-1)*dτ for i ∈ 1:Nτ]
println("τₑ = ", τ[end])
fun(τ,T) = dLdτ(τ,T, t, A, tᶜ,u, dt)

# Solve optimization problem
Tlb = 0*Tcut; Tub = 35
Tₒₚₜ = RK4_nln_sys(τ,T,fun,Tlb,Tub,ϵ)

# Optimal values
Tₒ = Tₒₚₜ[1]
for n ∈ 1:N
    if Tₒ[n] < Tcut
        A[n] = 0 
    end
end
Aₒ = A[:]
tᶜₒ = tᶜ[:]
lenT = length(Tₒ)

#############################################################################################
# Resulting Gaussian approximation of the envelope
G = gauss_fun(t, Aₒ, tᶜₒ, Tₒ)

# Propagation of resulting WGs
θ = angle.(𝓗)       # θ = ϕ-ωt
eθ = exp.(1im*θ)     # exp(iθ) = exp(-iωt)*exp(iϕ)
ηₜ = G.*eθ

## Instantaneous frequency
omH = zeros(Float64, M)
θᵤ = unwrap(θ)
for i ∈ 2:M-1
    omH[i] = (θᵤ[i+1] - θᵤ[i-1]) / (t[i+1]-t[i-1])
end

ωₒ = 2π ./ Tₒ
ωᵢ = zeros(Float64, lenT)
βᵢ = zeros(Float64, lenT)

for n ∈ 1:lenT
    eϕᵢ = eθ .* exp.(1im*ωₒ[n]*t)           # exp(iϕ) = exp(iθ)/exp(-iωt) = exp(iθ)*exp(iωt)
    ϕᵢ = real(-1im*log.(eϕᵢ))
    ϕᵢᵢ = unwrap(ϕᵢ)
    # ωᵢ[n] = (ϕᵢᵢ[end]-ϕᵢᵢ[1])/(t[end]-t[1])
    fit = linear_fit(t,ϕᵢᵢ)
    ωᵢ[n] = fit[2]    # ϕ = ω⁺*t + β
    βᵢ[n] = fit[1]
end

ω̃ = (ωᵢ.-ωₒ)[1]

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
eθ̃ = exp.(1im*(ω⁺*t .+ β)) .* exp.(1im*ωt)
# eθ̃ = exp.(1im*(ω⁺*t .+ β)) .* eωt
# ηₜ = G.*real(eθ̃)

#############################################################################################
# Normalization
Ω = ω̃/ωₚⱼ;      ω̅⁺ = ω⁺/ωₚⱼ;      ω̅ᵢ = ωᵢ./ωₚⱼ
G̅ = G ./ maximum(u)
t̅ = t ./ Tₚⱼ
T̅z = Tz / Tₚⱼ    
T̅IE = TIE ./ Tₚⱼ;   T̅DE = TDE ./ Tₚⱼ
T̅ₒ = Tₒ ./ Tₚⱼ;     A̅ₒ = Aₒ ./ Hₛⱼ;     t̅ᶜ = tᶜₒ ./ Tₚⱼ

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
Nₛ = Int(nextpow(2,M)/2)+1
SP_TOT = zeros(Float64, Nₛ) # Total spectrum of propagated Gaussian EWGs

pltEWG = plot(palette=:darkrainbow)
pltSPEC = plot(palette=:darkrainbow)
pltCONV = plot(palette=:darkrainbow)

for n ∈ 1:lenT
    ξ = -β .+ [π/2; π; π; π/2] # .+ zeros(Float64, lenT)
    # ξ = ωᵢ[n] * t .+ β
     
    # Gaussian EWG envelope
    gₙ = elem_wg(t̅*Tₚⱼ, A̅ₒ*Hₛⱼ, t̅ᶜ*Tₚⱼ, T̅ₒ*Tₚⱼ, n)  # Non-dimensional
    # gₙ = elem_wg(t, Aₒ, tₘ, Tₒ, n)
    # gₙ = elem_wg(t, Aₒ, tᶜₒ, Tₒ, n)
    # gₙ = elem_wg(t, Aₒ, zeros(Float64, lenT), Tₒ, n)

    # Propagated EWG
    WGT,_ = wave_group(1,T̅ₒ[n]*Tₚⱼ,1,2π/(Ω*ωₚⱼ),d,t̅*Tₚⱼ,0)
    # WGT,_ = wave_group(1,T̅ₒ[n]*Tₚⱼ,1,2π/(Ω*ωₚⱼ),d,t̅*Tₚⱼ,0) 
    # WGT,_ = wave_group(1,T̅ₒ[n]*Tₚⱼ,1,2π/ωₚⱼ,d,t̅*Tₚⱼ,0) 
    # WGT,_ = stokes_sum(1,Tₒ[n],1,2π/ω̃,d,t,0)
    ηₙ = gₙ .* WGT 
    # ηₙ = gₙ .* exp.(-1im * Ω*ωₚⱼ * t̅*Tₚⱼ) #.* exp.(1im * β)
    # ηₙ = gₙ .* exp.(-1im * ω̃ * t) #.* exp.(1im * β)
    # ηₙ = gₙ .* exp.(-1im * ωₒ[n] * t) .* exp.(1im * β)
    # ηₙ = gₙ .* exp.(1im*θ)
    FR, MAG, ang, df,_ = one_side_asp(real(ηₙ),t̅*Tₚⱼ) 

    # Propagated EWG - Focused WG
    # ηⁱₙ = gₙ .* exp.(-1im * Ω*ωₚⱼ * t̅*Tₚⱼ) 
    # ηⁱₙ = gₙ .* exp.(-1im * ω̃ * t) 
    # ηⁱₙ = gₙ .* exp.(-1im * ωₒ[n] * t)
    # ηⁱₙ = gₙ .* exp.(1im*θ)
    # ηₙ, MAG, FR = fcsd_wg(ηⁱₙ, t̅*Tₚⱼ, A̅ₒ*Hₛⱼ, t̅ᶜ*Tₚⱼ, βᵢ[n], n)
    # ηₙ, MAG, FR = fcsd_wg(ηⁱₙ, t̅*Tₚⱼ, A̅ₒ*Hₛⱼ, t̅ᶜ*Tₚⱼ, ξ[n], n)        

    TOT[:] = TOT[:] .+ ηₙ
    SP_TOT[:] = SP_TOT[:] .+ MAG

    plot!(pltEWG, t̅*Tₚⱼ, real(ηₙ), xlab = "t [s]", ylab = "A [m]", lab = "g$(n)", lw=2, legend=:outerbottom, legendcolumns=lenT)
    plot!(pltSPEC, FR, MAG, xlab = L"f [Hz]", ylab = L"A [m^2 s]", lab = "g$(n)", lw=2, legend=:outerbottom, legendcolumns=lenT)
    plot!(pltCONV, Tₒₚₜ[2][:,n], xlab = "iteration",  ylab = L"T_o^i [s]", xscale=:log10, lab = "To[$(n)]", lw=2, legend=:outerbottom, legendcolumns=lenT)

    for i ∈ 1:M
        if abs(ηₙ[i]) < 1e-6
            ηₙ[i] = 0
        end
    end                             

    fid_EWG = "EWG_$n" # File name
    open(fid_res*fid_EWG, "w")
    # writedlm(fid_res*fid_EWG, [round.(t*1e6)./1e6 round.(gₙ*1e6)./1e6 round.(real(ηₙ*1e6))./1e6 round.(angle.(ηₙ)*1e6)./1e6], '\t')
    writedlm(fid_res*fid_EWG, [round.(t*1e6)./1e6 round.(gₙ*1e6)./1e6 round.(WGT*1e6)./1e6 round.(ηₙ*1e6)./1e6], '\t')
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
plot!(t̅*Tₚⱼ, real(TOT), lab = L"\sum EWGs", lw=2, line=:dot)

plt_eta2 = plot(xlab=L"t~[s]", ylab=L"\eta ~[m]", palette=[cb[4];cb[8]]) 
plot!(t, real(ηₜ), lab=L"\eta_G(t)", lw=2, line=:dot)
plot!(t̅*Tₚⱼ, real(TOT), lab = L"\sum EWGs", lw=2)

# Spectra of resulting signals
plt_spec = plot(palette=[cb[11];cb[4];cb[8];cb[1]])
plot!(fr_η, mag_η, xlab = "f [Hz]", ylab = L"A [m^2 s]", lab = "Original", lw=2)
freq, mag, _ = one_side_asp(real(ηₜ), t)
plot!(freq, mag, xlab = "f [Hz]", ylab = L"A [m^2 s]", lab = L"\eta_{new}", lw=2)
freq, mag, _ = one_side_asp(real(TOT), t̅*Tₚⱼ)
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

# Display plots
display(plt_sigenv)
display(plt_env)
display(plt_G)
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
plot!(tᵢ,etaG, lw=2, lab=L"G (t)")

EWGsum = zeros(Float64, Lᵢ)
EWGsum[istart:iend] = real(TOT)[:]
plot!(tᵢ,EWGsum, lw=2, lab=L"\sum g_i (t)")

ηTOT = fcsd_wg(freq, SP_TOT, t, Aₒ[1], 0, 0)
etaFCSD = zeros(Float64, Lᵢ)
etaFCSD[istart:iend] = ηTOT[:]
display(plot!(tᵢ, etaFCSD, lw=2, lab="Focused wave"))

# Save figures
savefig(plt_G, fid_res*"env_Ginterp.svg");  savefig(plt_G, fid_res*"env_Ginterp.png")
savefig(pltEWG, fid_res*"EWGs.svg");    savefig(pltEWG, fid_res*"EWGs.png")
savefig(plt_eta0, fid_res*"new_sig_comp.svg");   savefig(plt_eta0, fid_res*"new_sig_comp.png")
savefig(plt_pars, fid_res*"EWGpars.svg");   savefig(plt_pars, fid_res*"EWGpars.png")

#############################################################################################
# Write results in text files
t = round.(t*1e6)./1e6
ηₜ = round.(ηₜ*1e6)./1e6
TOT = round.(TOT*1e6)./1e6
tᶜₒ = round.(tᶜₒ*1e6)./1e6
Aₒ = round.(Aₒ*1e6)./1e6
Tₒ = round.(Tₒ*1e6)./1e6

open(fid_res*"etaG", "w")
# writedlm(fid_res*"etaG", [t real(ηₜ) angle.(ηₜ)], '\t')
writedlm(fid_res*"etaG", [tᵢ etaG], '\t')

open(fid_res*"EWGsum", "w")
# writedlm(fid_res*"EWGsum", [t real(TOT) angle.(TOT)], '\t')
writedlm(fid_res*"EWGsum", [tᵢ EWGsum], '\t')

open(fid_res*"etaFCSD", "w")
writedlm(fid_res*"etaFCSD", [tᵢ etaFCSD], '\t')

head = ["tᶜₒ [s]" "Tₒ [s]" "Aₒ [m]"]
open(fid_res*"EWGpars", "w")
writedlm(fid_res*"EWGpars", [head; tᶜₒ Tₒ Aₒ], '\t')    

# head = ["t̅ᶜ [-]" "T̅ₒ [-]" "A̅ₒ [-]" "βᵢ [rad]" "ω̅ᵢ [-]" "T̅DE [-]" "t̅ₘ"]
head = ["t̅ᶜ [-]" "T̅ₒ [-]" "A̅ₒ [-]" "βᵢ [rad]" "ω̅ᵢ [-]" "T̅ₛ" "T̅ₕ" "T̅ₚ"]
open(fid_res*"EWG_norm_pars", "w")
# row = round.([t̅ᶜ T̅ₒ A̅ₒ βᵢ ω̅ᵢ T̅DE t̅ₘ]*1e6)./1e6
row = round.([t̅ᶜ T̅ₒ A̅ₒ βᵢ ω̅ᵢ Tₒ./T₀ₛ² Tₒ./T₀ₕ Tₒ./T₀ₚ]*1e6)./1e6
writedlm(fid_res*"EWG_norm_pars", [head; row], '\t')

head = ["Hₛ [m]" "Tₚ [s]" "$std_fac*σₙ" "ηₘₐₓ [m]" "uₘₐₓ [m]" "t̅ₑᵥ [-]" "Ω [-]" "β [rad]" "ω̅⁺ [-]" "T̅z [-]"]
open(fid_res*"G_pars", "w")
row = round.([Hₛⱼ Tₚⱼ MinPeakVal ηHP maximum(u) t̅[end]-t̅[1] Ω β ω̅⁺ T̅z]*1e6)./1e6
writedlm(fid_res*"G_pars", [head; row], '\t')

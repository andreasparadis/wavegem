# Load necessary modules
using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles
using FFTW, DSP, LsqFit, Statistics

# Include necessary scripts for functions
include("signal_processing.jl")
include("peak_detect.jl")
include("runge_kutta.jl")
include("Gauss_Reg_t.jl")
include("directories.jl")

g ::Float64 = 9.81
# Read signal from text file in library folder
dir0::String = "SE/LAF/2/"
fid0::String = dir0*"Decomposition/events/"
fid_ev::String = "event_1" # File name
fid = fid0*fid_ev # Path to file
cont = Array{Float64}(undef,0,2)

open("library/"*fid, "r")
cont = readdlm("library/"*fid, '\t', Float64, skipstart=1, '\n')

t = round.((cont[:,1] .- cont[1,1])*100)/100 # x-coordinates vector
η = cont[:,2] # Surface elevation
M = length(t);      Øᴹ = zeros(Float64,M) 
dt = (t[end]-t[1])/(M-1)

# Hilbert transform and envelope of surface elevation
𝓗 = hilbert(η);   u = abs.(𝓗)

# Peak detection based on 2nd derivative of signal
# tᶜᵢ, Aᵢ, uₜ, uₜₜ, i⁺ = maxima(u, t)

fcut = 0.625 # Cut-off frequency [Hz]
Tcut = 1/fcut
MinPeakVal, MinPeakDist = 1*std(η), Tcut
Aᵢ, tᶜᵢ, i⁺, uₜ, uₜₜ = peaks_max_ext(u,t, MinPeakVal, MinPeakDist)

# Additional filter on peaks
# 1: None , 2: 2nd derivative
Filter = 2
rng = Int(round(Tcut/4/dt)+1)
len = length(tᶜᵢ)
tᶜ = [tᶜᵢ[1]]
A = [Aᵢ[1]]

if Filter == 1
    tᶜ = tᶜᵢ[:]
    A = Aᵢ[:]
elseif Filter == 2
    for i ∈ 2:len
        if sign(uₜ[i⁺[i]-rng]) != sign(uₜ[i⁺[i]+rng])
            push!(tᶜ,tᶜᵢ[i])
            push!(A,Aᵢ[i])
        end
    end
end

N = length(tᶜ);     Øᴺ = zeros(Float64,N)

# Perform the least squares fitting
## Initial condition for L based on ∂ₜ²g(tᶜ) = ∂ₜ²u(tᶜ)
T = Øᴺ

for n ∈ 1:N
    i = i⁺[n]
    T[n] = sqrt.(-2*A[n]./uₜₜ[i])
end

# Definition of the ODE's dLₘ  
τ = collect(range(0,1,M))

fun(τ,T) = dLdτ(τ,T, t, A, tᶜ,u, dt)
Tₒₚₜ = RK4_nln_sys(τ,T,fun)
Tₒ = Tₒₚₜ[1]
Aₒ = A
tᶜₒ = tᶜ

# Low-pass filter of optimal values
# k_min = 0.5 * 2π/t[end];    T_max = 2π / sqrt(g*k_min)
# Tₒ = [0.0];     tᶜₒ = [0.0];    Aₒ = [0.0]

# for i ∈ 1:N
#     if Tₒₚₜ[1][i] > Tcut # && Tₒₚₜ[1][i] < T_max
#         push!(Tₒ, Tₒₚₜ[1][i])
#         push!(tᶜₒ, tᶜ[i])
#         push!(Aₒ, A[i])
#     end
# end

# Tₒ = Tₒ[2:end]; 
# tᶜₒ = tᶜₒ[2:end]
# Aₒ = Aₒ[2:end]

# Sort in descending order the pairs (Aₒ,Tₒ) based on Aₒ
lenT = length(Tₒ)
for i in 1:lenT
    for j in 1:lenT-i
        if Aₒ[j] < Aₒ[j+1]
            Tₒ[j], Tₒ[j+1] = Tₒ[j+1], Tₒ[j]
            Aₒ[j], Aₒ[j+1] = Aₒ[j+1], Aₒ[j]
            tᶜₒ[j], tᶜₒ[j+1] = tᶜₒ[j+1], tᶜₒ[j]
        end
    end
end

G = gauss_fun(t, Aₒ, tᶜₒ, Tₒ)

# Signal and envelope
lbl_x = L"t~[s]";   lbl_y = L"[m]";   llbls = [L"η(t)" L"u(t)"]
plt = plot(t, [η u], xlab=lbl_x, ylab=lbl_y, lab=llbls, lw=[1 2])
display(plt)

# Scatter plot of extrema
plot(t, [u G], xlab = L"t~[s]", ylab = L"[m]", lab = [L"u(t)" L"G(t)"], lw=[2 1])
display(plot!(tᶜ, A, seriestype=:scatter, ms=2, mc=:red, lab = "filtered peaks"))
# display(plot!(tᶜₒ, Aₒ, seriestype=:scatter, ms=2, mc=:green, lab = "final filtered peaks"))
savefig("library/"*fid0*"env_Ginterp.svg")
savefig("library/"*fid0*"env_Ginterp.png")

# Propagation of resulting WGs
ηₜ = real(G.*exp.(-1im*angle.(𝓗)))


lbl_x = L"t~[s]";   lbl_y = L"\eta ~[m]";   llbls = L"\eta_{new}(t)" 
plt = plot(t, ηₜ, xlab=lbl_x, ylab=lbl_y, lab=llbls, lw=2)
plot!(t,η, lab = "Original signal")
display(plt)
savefig("library/"*fid0*"new_sig_comp.svg")
savefig("library/"*fid0*"new_sig_comp.png")

open("library/"*fid0*"etaG", "w")
writedlm("library/"*fid0*"etaG", [t ηₜ], '\t')

# Plot the group parameters in ascending order

display(plot(Tₒ, Aₒ, seriestype=:scatter, ms=2, mc=:red, xlab = L"T_o~[s]", ylab = L"A_o~[m]", lab = "Group parameters (sorted)"))
savefig("library/"*fid0*"EWGpars.svg")
savefig("library/"*fid0*"EWGpars.png")

fid_pars = "EWGpars" # File name
head = ["tᶜₒ [s]" "Tₒ [s]" "Aₒ [m]"]
open("library/"*fid0*fid_pars, "w")
writedlm("library/"*fid0*fid_pars, [head; tᶜₒ Tₒ Aₒ], '\t')

# lenEWG = 1000
# EWGs = zeros(Float64, M)
TOT = zeros(Float64, M)
pltEWG = plot()

for n ∈ 1:lenT
    # tEWG = range(-Tₒ[n]/2, Tₒ[n]/2, lenEWG)
    # EWGs = real(Aₒ[n].*exp.(-1im* 2*pi/Tₒ[n] .*tEWG))
    EWGs = elem_wg(t, Aₒ, tᶜₒ, Tₒ, n) .* real(exp.(-1im* 2*pi/Tₒ[n] * t))

    TOT[:] = TOT[:] + EWGs
    plot!(t, EWGs, xlab = "t [s]", ylab = "A [m]", lab = "g$(n)", legend=:outerbottom, legendcolumns=7)

    for i ∈ 1:length(EWGs)
        if abs(EWGs[i]) < 1e-6
            EWGs[i] = 0
        end
    end

    fid_EWG = "EWG_$n" # File name
    open("library/"*fid0*fid_EWG, "w")
    writedlm("library/"*fid0*fid_EWG, [t EWGs], '\t')
end

display(pltEWG)
savefig("library/"*fid0*"EWGs.svg")
savefig("library/"*fid0*"EWGs.png")

plot(plt)
display(plot!(t,TOT, lab = "sum of EWGs"))
savefig("library/"*fid0*"EWGsum.svg")
savefig("library/"*fid0*"EWGsum.png")

fid_GS = "EWGsum" # File name
open("library/"*fid0*fid_GS, "w")
writedlm("library/"*fid0*fid_GS, [t TOT], '\t')
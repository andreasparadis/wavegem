# Load necessary modules
using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles
using FFTW, DSP, LsqFit

# Include necessary scripts for functions
include("signal_processing.jl")
include("peak_detect.jl")
include("runge_kutta.jl")
include("Gauss_Reg.jl")

g ::Float64 = 9.81
# Read signal from text file in library folder
fid::String = "rand_seas/JS1_H3_T8_eta-t.txt" # Path to file
cont = Array{Float64}(undef,0,2)

open("library/"*fid, "r")
cont = readdlm("library/"*fid, '\t', Float64, skipstart=1, '\n')

x = cont[:,1] # x-coordinates vector
η = cont[:,2] # Surface elevation
M = length(x);      Øᴹ = zeros(Float64,M) 

# Hilbert transform and envelope of surface elevation
𝓗 = hilbert(η);   u = abs.(𝓗)

# Peak detection based on 2nd derivative of signal
xᶜᵢ, Aᵢ, uₓ, uₓₓ, i⁺ = maxima(u, x)

dt = (x[end]-x[1])/(length(x)-1)
fcut = 2 # Cut-off frequency [Hz]
Tcut = 1/fcut
rng = Tcut/dt
len = length(xᶜᵢ)
xᶜ = [xᶜᵢ[1]]
A = [Aᵢ[1]]

i = 1
while i < len
    if xᶜᵢ[i+1]-xᶜᵢ[i] < Tcut
        # push!(xᶜ,(xᶜᵢ[i-1]+xᶜᵢ[i+1])/2)
        # push!(A,(Aᵢ[i-1]+Aᵢ[i+1])/2)
        global i += 2
    else
        push!(xᶜ,xᶜᵢ[i+1])
        push!(A,Aᵢ[i+1])
        global i += 1
    end
end

N = length(xᶜ);     Øᴺ = zeros(Float64,N)

# Gaussian Regression
# g, G, dGdL, If₁ = (Øᴹ[:] for _ = 1:4)

# Perform the least squares fitting
## Initial condition for L based on ∂ₓ²g(xᶜ) = ∂ₓ²u(xᶜ)
L = Øᴺ

for n ∈ 1:N
    i = i⁺[n]
    L[n] = sqrt.(-2*A[n]./uₓₓ[i])
end

# Definition of the ODE's dLₘ  
τ = collect(range(0,1,M))
dx = (x[end]-x[1])/(M-1)

fun(τ,L) = dLdτ(τ,L)
Lₒₚₜ = RK4_nln_sys(τ,L,fun)
G = gauss_fun(x, A, xᶜ, Lₒₚₜ[1])

# # Low-pass filter
# fcut = 5 # Cut-off frequency [Hz]
# Tcut = 1/fcut
# k_min = 0.5 * 2π/x[end]
# T_max = 2π / sqrt(g*k_min)
# Lₒ = [0.0]
# xᶜₒ = [0.0]
# Aₒ = [0.0]

# for i ∈ 1:N
#     if Lₒₚₜ[1][i] > Tcut && Lₒₚₜ[1][i] < T_max
#         push!(Lₒ, Lₒₚₜ[1][i])
#         push!(xᶜₒ, xᶜ[i])
#         push!(Aₒ, A[i])
#     end
# end

# Lₒ = Lₒ[2:end]
# xᶜₒ = xᶜₒ[2:end]
# Aₒ = Aₒ[2:end]
# G = gauss_fun(x, Aₒ, xᶜₒ, Lₒ)

# # Group 1
# L1 = [ Lₒ[18]; L[19] ]
# x1 = [xᶜₒ[18] xᶜₒ[19]]
# A1 = [ Aₒ[18] Aₒ[19] ]

# ηₜ = zeros(Float64,M)
# for i ∈ 1:2
#     for m ∈ 1:M
#         ηₜ[m] = ηₜ[m] + A1[i] * cos(2π/L1[i]*(x[m]-x1[i]))
#     end
# end

# Signal and envelope
lbl_x = L"x~[m]";   lbl_y = L"[m]";   llbls = [L"η(x)" L"u(x)"]
plt = plot(x, [η u], xlab=lbl_x, ylab=lbl_y, lab=llbls, lw=[1 2])
display(plt)

# Scatter plot of extrema
plot(x, [u G], xlab = L"x~[m]", ylab = L"[m]", lab = [L"u(x)" L"G(x)"], lw=[2 1])
display(plot!(xᶜ, A, seriestype=:scatter, ms=2, mc=:red, lab = "maxima"))
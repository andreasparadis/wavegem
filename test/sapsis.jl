# Load necessary modules
using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles
using FFTW, DSP, LsqFit

# Include necessary scripts for functions
include("signal_processing.jl")
include("peak_detect.jl")
include("runge_kutta.jl")

# Read signal from text file in library folder
fid::String = "rand_seas/JS1_H3_T8_eta-x.txt" # Path to file
cont = Array{Float64}(undef,0,2)

open("library/"*fid, "r")
cont = readdlm("library/"*fid, '\t', Float64, skipstart=1, '\n')

x = cont[:,1] # x-coordinates vector
η = cont[:,2] # Surface elevation

# Hilbert transform and envelope of surface elevation
𝓗 = hilbert(η);   u = abs.(𝓗)

# Peak detection based on 2nd derivative of signal
xᶜ, A, uₓ, uₓₓ, i⁺ = maxima(u, x)

M = length(x);      Øᴹ = zeros(Float64,M) 
N = length(xᶜ);     Øᴺ = zeros(Float64,N)

# Gaussian Regression
# g, G, dGdL, If₁ = (Øᴹ[:] for _ = 1:4)

## Define the Elementary Gaussian model function
function elem_wg(x, a, xᶜ, L, n)
    g = a[n] * exp.(-((x.-xᶜ[n])/L[n]).^2)

    return g
end

# Define the Gaussian model function
function gauss_fun(x, a, xᶜ, L)
    N = length(xᶜ)
    G = sum(elem_wg(x, a, xᶜ, L, n) for n in 1:N)

    return G
end

# Define the objective function for least squares fitting
function obj_fun(L, a, xᶜ, x, u)
    M = length(x)
    dx = (x[end]-x[1]/M)
    J = 1/2 * sum((gauss_fun(x, a, xᶜ, L) - u).^2 *dx for m ∈ 1:M) 
    return J
end

# Perform the least squares fitting

## Initial condition for L based on ∂ₓ²g(xᶜ) = ∂ₓ²u(xᶜ)
L = Øᴺ

for n ∈ 1:N
    i = i⁺[n]
    L[n] = sqrt.(-2*A[n]./uₓₓ[i])
end

## Apply fit
# J(L) = obj_fun(L, A, xᶜ, x, u)

fit_res = curve_fit((ptest, xtest) -> obj_fun(Øᴺ, A, xᶜ, x, u), x, u, L)
Lₒₚₜ = fit_res.param

G = gauss_fun(x, A, xᶜ, Lₒₚₜ)

# Signal and envelope
lbl_x = L"x~[m]";   lbl_y = L"[m]";   llbls = [L"η(x)" L"u(x)"]
plt = plot(x, [η u], xlab=lbl_x, ylab=lbl_y, lab=llbls, lw=[1 2])
display(plt)

# Scatter plot of extrema
plot(x, [u G], xlab = L"x~[m]", ylab = L"[m]", lab = [L"u(x)" L"G(x)"], lw=[2 1])
display(plot!(xᶜ, A, seriestype=:scatter, ms=2, mc=:red, lab = "maxima"))
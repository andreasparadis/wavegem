# Load necessary modules
using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles
using FFTW, DSP, CurveFit, Distributions

# Include necessary scripts for functions
include("signal_processing.jl")
include("peak_detect.jl")

N, M::Int32 = 100, 1000

x = range(0,200,M)
Lₓ = rand(Uniform(2,50),N)
A = rand(Uniform(0,5),N)
xᶜ = rand(Uniform(0,200),N)

g = zeros(Float64,N)   
G = zeros(Float64,M)

for m ∈ 1:M
    for n ∈ 1:N
        g[n] = A[n] * exp.(-((x[m]-xᶜ[n])/Lₓ[n]).^2)
        G[m] = G[m] + g[n]
    end
end

plot(x, G)


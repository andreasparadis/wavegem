# Plot the final surface elevation from the OW3D output text file "eta0_coeffs"
using Plots
using DataFrames
using LaTeXStrings
using LinearAlgebra
using DelimitedFiles
using FFTW, FourierAnalysis, DSP

include("signal_processing.jl")
include("jonswap.jl")

fid::String = "run_8/eta0_coeffs.txt"
cont = Array{Float64}(undef,0)
Hₛ, Tₚ, Tᵢ, Tₑ, γ::Float64 = 3.0, 8.0, 0.4, 1000, 3.3 # Real scale
# Hₛ, Tₚ, Tᵢ, Tₑ, γ::Float64 = 0.0375, 0.8944, 0.2, 3, 3.3 # Experimental scale

fⱼ,Sⱼ = spectrum(Hₛ,Tₚ,γ,Tᵢ,Tₑ,1000) # Generate JONSWAP spectrum for pair (Hₛ,Tₚ) ∈ [Tᵢ,Tₑ]

open("library/"*fid, "r")
cont = readdlm("library/"*fid, '\t', Float64, skipstart=1, '\n')

N = length(cont[:,1])
fₙ = zeros(Float64,N,1)
cₙ = zeros(Complex{Float64},N,1)
Cₙ = zeros(Float64,N,1)
Re_cₙ = zeros(Float64,N,1)
Im_cₙ = zeros(Float64,N,1)

fₙ = cont[:,1]
Cₙ = cont[:,2]
Re_cₙ = cont[:,3]
Im_cₙ = cont[:,4]

ωₙ = 2π * fₙ

dt = 0.05
tₑ = dt*N
t = range(0, tₑ, N)

η = zeros(Complex{Float64},N,1)

for i ∈ 1:N
    cₙ[i] = complex(Re_cₙ[i], Im_cₙ[i])
    for j ∈ 1:N
        η[j] = η[j] + cₙ[i] * exp(ωₙ[i]*t[j]*im)
    end
end

plot(t, real(η), xlab = L"t~[s]", ylab = L"\eta~[]", lab = L"Surface~elevation", lw=1)

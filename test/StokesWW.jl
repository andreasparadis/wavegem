using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles, FFTW
using Statistics, Random, Distributions
using Dates

include("wave_theory.jl")
include("signal_processing.jl")
include("jonswap.jl")
include("peak_detect.jl")

H₁, T₁::Float64 = 10.0, 8.0
H₂, T₂::Float64 = 6.0, 8.0
d, ρ, g::Float64 = 10.0, 1025.0, 9.81

ω₁ = 2π / T₁ 
λ₁⁰ = g/(2*π) * T₁^2
print("$λ₁⁰ \n")
λ₁,_ = dispersion(T₁,d)
υ₁ = λ₁/T₁;     κ₁ = 2π/λ₁

ω₂ = 2π / T₂ 
λ₂⁰ = g/(2*π) * T₂^2
print("$λ₂⁰ \n")
λ₂,_ = dispersion(T₂,d)   
υ₂ = λ₂/T₂;     κ₂ = 2π/λ₂

Nₜ = 1000
tₑ = 2 * max(T₁, T₂)
dt = tₑ/(Nₜ-1)
t = zeros(Float64,Nₜ)
[t[i] = (i-1)*dt for i ∈ 1:Nₜ]

η̂₁ = H₁/2
η̂₂ = H₂/2
θ₁ = -ω₁ * t
θ₂ = -ω₂ * t
η₁ = zeros(Float64,Nₜ)
η₂ = zeros(Float64,Nₜ)

η₁ = η₁ .+ η̂₁ * cos.(θ₁) .+ η̂₁^2 * κ₁/4 * cosh(κ₁*d)*(2+cosh(2κ₁*d))/sinh(κ₁*d)^3 * cos.(2*θ₁) 
η₂ = η₂ .+ η̂₂ * cos.(θ₂) .+ η̂₂^2 * κ₂/4 * cosh(κ₂*d)*(2+cosh(2κ₂*d))/sinh(κ₂*d)^3 * cos.(2*θ₂) 

Nₓ = 1000
x = collect(range(0,2*λ₁,Nₓ));
ϕ₁ = κ₁*x
η₁ˣ = zeros(Float64,Nₓ)

η₁ˣ = η₁ˣ .+ η̂₁ * cos.(ϕ₁) .+ η̂₁^2 * κ₁/4 * cosh(κ₁*d)*(2+cosh(2κ₁*d))/sinh(κ₁*d)^3 * cos.(2*ϕ₁) 

display(plot(t, [η₁ η̂₁*cos.(θ₁)], lab=["Stokes" "Linear"]))
display(plot(t, [η₂ η̂₂*cos.(θ₂)], lab=["Stokes" "Linear"]))
plot(x, η₁ˣ, aspect_ratio=1)

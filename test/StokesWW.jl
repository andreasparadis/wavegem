using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles, FFTW
using Statistics, Random, Distributions
using Dates

include("wave_theory.jl")
include("signal_processing.jl")
include("jonswap.jl")
include("peak_detect.jl")

H₁, T₁::Float64 = 5.139948, 8.523
H₂, T₂::Float64 = 5.139948, 5.081
# H₁, T₁::Float64 = 1, 1
# H₂, T₂::Float64 = 1, 1/sqrt(2)
d, ρ, g::Float64 = 100.0, 1025.0, 9.81

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

λw = 4π/(κ₁+κ₂)
Tw = 4π/(ω₁+ω₂)
λg = 4π/(κ₁-κ₂)
Tg = 4π/(ω₁-ω₂)

Nₜ = 1000
tₑ = 2*round(Tg)
dt = tₑ/(Nₜ-1)
t = zeros(Float64,Nₜ)
[t[i] = (i-1)*dt-Tg for i ∈ 1:Nₜ]

η̂₁ = H₁/2
η̂₂ = H₂/2
θ₁ = -ω₁ * t
θ₂ = -ω₂ * t

α₁ = coth(κ₁*d)
α₂ = coth(κ₂*d)
Cnum = (2*ω₁*ω₂*(ω₁-ω₂)*(1+α₁*α₂) + ω₁^3*(α₁^2-1) - ω₂^3*(α₂^2-1))*(ω₁-ω₂)*(α₁*α₂-1)
Cden = ω₁^2*(α₁^2-1) - 2*ω₁*ω₂*(α₁*α₂-1) + ω₂^2*(α₂^2-1) + (ω₁^2+ω₂^2) - ω₁*ω₂*(α₁*α₂+1) 
C = Cnum/Cden
Dnum = (2*ω₁*ω₂*(ω₁+ω₂)*(α₁*α₂-1) + ω₁^3*(α₁^2-1) + ω₂^3*(α₂^2-1))*(ω₁+ω₂)*(α₁*α₂+1)
Dden = ω₁^2*(α₁^2-1) - 2*ω₁*ω₂*(α₁*α₂+1) + ω₂^2*(α₂^2-1) - (ω₁^2+ω₂^2) + ω₁*ω₂*(α₁*α₂-1) 
D = Dnum/Dden

η₁ = zeros(Float64,Nₜ)
η₂ = zeros(Float64,Nₜ)

η₁ = η₁ .+ η̂₁ * cos.(θ₁) .+ η̂₁^2 * κ₁/4 * cosh(κ₁*d)*(2+cosh(2κ₁*d))/sinh(κ₁*d)^3 * cos.(2*θ₁) 
η₂ = η₂ .+ η̂₂ * cos.(θ₂) .+ η̂₂^2 * κ₂/4 * cosh(κ₂*d)*(2+cosh(2κ₂*d))/sinh(κ₂*d)^3 * cos.(2*θ₂) 
η = η₁ .+ η₂ .+ α₁*α₂/(2*g) * (C*cos.(θ₁-θ₂) - D*cos.(θ₁+θ₂))
ηlin = (η̂₁+η̂₂) * cos.(2π/Tg *t).*cos.(2π/Tw *t)

Nₓ = 1000
x = collect(range(0,2*λ₁,Nₓ));
ϕ₁ = κ₁*x
η₁ˣ = zeros(Float64,Nₓ)

η₁ˣ = η₁ˣ .+ η̂₁ * cos.(ϕ₁) .+ η̂₁^2 * κ₁/4 * cosh(κ₁*d)*(2+cosh(2κ₁*d))/sinh(κ₁*d)^3 * cos.(2*ϕ₁) 

display(plot(t, [η₁ η̂₁*cos.(θ₁)], lab=["Stokes" "Linear"]))
display(plot(t, [η₂ η̂₂*cos.(θ₂)], lab=["Stokes" "Linear"]))
display(plot(x, η, aspect_ratio=1))
display(plot(t, [η ηlin]))
display(plot(t, [ηlin (η̂₁+η̂₂)/2*cos.(2π/Tg *t) (η̂₁+η̂₂)/2*cos.(2π/Tw *t)]))

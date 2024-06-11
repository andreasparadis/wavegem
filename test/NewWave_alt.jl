# NewWave Design Wave Group
using Plots
using DataFrames
using LaTeXStrings
using LinearAlgebra
using DelimitedFiles
using FFTW
using Statistics
using Random
using Distributions

include("signal_processing.jl")
include("jonswap.jl")
include("wave_theory.jl")

Hₛ, Tₚ, Tᵢ, Tₑ, γ, g::Float64 = 5.3, 8.0, 1.6, 25, 3.3, 9.81
x₀::Float64 = 0.0

dt = 0.1;   tend = 100
M = Int64(round(tend/dt))+1
t = zeros(Float64,M)
t₀ = tend/2

[t[i] = (i-1)*dt for i ∈ 1:M]

ω⁻ = 2π/Tₑ
ω⁺ = 2π/Tᵢ
dω = 2π/((M-1)*dt)
Nₛ = Int64(round((ω⁺-ω⁻)/dω))

fⱼ,Sⱼ = spectrum(Hₛ,Tₚ,γ,Tᵢ,Tₑ, Nₛ) # Generate JONSWAP spectrum for pair (Hₛ,Tₚ) ∈ [Tᵢ,Tₑ]

ωₚ = 2π / Tₚ
kₚ = ωₚ^2/g
λₚ = 2π / kₚ
ω = 2π * fⱼ
k = 1/g * ω.^2 
df = abs(fⱼ[end]-fⱼ[1])/(Nₛ-1)
m₀ = Hₛ^2/16

x = range(0,20*λₚ,M)
x₀ = x[end]/2
χ = x .- x₀

ηₓ = zeros(Float64,M)
ηₜ = zeros(Float64,M)
η̂ = sqrt.(2*Sⱼ*df)
# α = 8.9*sqrt(m₀)
α = sqrt(2*m₀*log(Nₛ))

# Spectrum energy through integration == Variance m₀
E = zeros(Float64,1)
for i ∈ 1:Nₛ
    E[1] = E[1] + Sⱼ[i]*df
end

for i ∈ 1:Nₛ
    for m ∈ 1:M
        ηₓ[m] = ηₓ[m] + Sⱼ[i]*df * cos(-ω[i]*(t[m] - t₀))
    end
end

ηₓ = α/m₀ * ηₓ

for i ∈ 1:Nₛ
    for m ∈ 1:M
        ηₜ[m] = ηₜ[m] + Sⱼ[i]*df * cos(k[i]*χ[m])
    end
end

ηₜ = α/m₀ * ηₜ

freq, mag, _ = one_side_asp(ηₓ,t)

display(plot(fⱼ, Sⱼ, xlab = L"f~[Hz]", ylab = L"\tilde{\eta}~[]", lab = L"JONSWAP", lw=3))
# savefig("JONSWAP.svg")
display(plot(fⱼ, η̂, xlab = L"f~[Hz]", ylab = L"\tilde{\eta}~[]", lab = L"Component Amplitudes", lw=3))
plot(t, ηₓ./α, xlab = L"\tau~[s]", ylab = L"\frac{\eta}{\alpha}~[-]", lw=3)
display(plot!(twiny(),χ, [ηₓ./α ηₜ./α], xlab = L"\chi~[m]", lab = [L"(\chi,\tau)=(0,\tau)" L"(\chi,\tau)=(\chi,0)"], ls=[:solid :dot], lw=3))
# savefig("NewWave.svg")
plot(freq, mag, xlab = L"f~[Hz]", ylab = L"Magnitude", lab = "Single-sided Spectrum", lw=3)

# fid0 = "library/SE/LAF/2/Decomposition/events/NewWave2"
# open(fid0, "w")
# writedlm(fid0, [t ηₓ], '\t')
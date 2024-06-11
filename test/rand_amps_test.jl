using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles, FFTW
using Statistics, Random, Distributions

include("wave_theory.jl")
include("signal_processing.jl")
include("jonswap.jl")

d, Hₛ, Tₚ, γ::Float64 = 100.0, 3, 8.0, 3.3
ρ, g::Float64 = 1025.0, 9.81

M = 1000
Tᵢ = 1.6
Tₑ = 25.0

fⱼ,Sⱼ = spectrum(Hₛ,Tₚ,γ,Tᵢ,Tₑ,M) # Generate JONSWAP spectrum for pair (Hₛ,Tₚ) ∈ [Tᵢ,Tₑ]
dfⱼ = abs(fⱼ[end]-fⱼ[1])/(M-1)
η̂ = sqrt.(2*Sⱼ*dfⱼ)


A = abs.(rand(Normal(),M))
cvar = cov(A, η̂)
sigA = std(A)
sigB = std(η̂)

# Pearson correlation coefficient
rho = cvar ./ (sigA * sigB)
print("ρ=", rho, "\n")

amp = rand(Normal(),M) + 1im*rand(Normal(),M) # Random complex vector (Normal distribution, 0 mean, 1 variance)
Aᴬ = abs.(amp)
cvar = cov(Aᴬ, η̂)
sigA = std(Aᴬ)
rho = cvar ./ (sigA * sigB)
print("ρ=", rho, "\n")

Random.seed!(Int64(round(rand(1)[1]*1000)));    distr = Rayleigh(1); # Rayleigh distribution with unit scale
Aᴿ = rand(distr,M)
cvar = cov(Aᴿ, η̂)
sigA = std(Aᴿ)
rho = cvar ./ (sigA * sigB)
print("ρ=", rho, "\n")

plt = plot(fⱼ, A.*η̂, xlab = L"f~[Hz]", ylab = L"Amplitude~[m]", lab = "|Normal|", lw=1, line=:dash)
plot!(fⱼ, Aᴬ.*η̂, xlab = L"f~[Hz]", ylab = L"Amplitude~[m]", lab = "|complex(Normal)|", lw=1, line=:dot)
plot!(fⱼ, Aᴿ.*η̂, xlab = L"f~[Hz]", ylab = L"Amplitude~[m]", lab = "Rayleigh", lw=1)
plot!(fⱼ, η̂, xlab = L"f~[Hz]", ylab = L"Amplitude~[m]", lab = "JONSWAP", lw=3)
display(plt)

plt = plot(fⱼ, A.*η̂, xlab = L"f~[Hz]", ylab = L"Amplitude~[m^2/Hz]", lab = "|Normal|", lw=1)
plot!(fⱼ, η̂, xlab = L"f~[Hz]", ylab = L"Amplitude~[m]", lab = "JONSWAP", lw=3)
display(plt)

plt = plot(fⱼ, Aᴬ.*η̂, xlab = L"f~[Hz]", ylab = L"Amplitude~[m^2/Hz]", lab = "|complex(Normal)|", lw=1)
plot!(fⱼ, η̂, xlab = L"f~[Hz]", ylab = L"Amplitude~[m]", lab = "JONSWAP", lw=3)
display(plt)

plt = plot(fⱼ, Aᴿ.*η̂, xlab = L"f~[Hz]", ylab = L"Amplitude~[m]", lab = "Rayleigh", lw=1)
plot!(fⱼ, η̂, xlab = L"f~[Hz]", ylab = L"Amplitude~[m]", lab = "JONSWAP", lw=3)
display(plt)
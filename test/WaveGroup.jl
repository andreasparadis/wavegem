using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles, FFTW
using Statistics, Random, Distributions

include("signal_processing.jl")
include("jonswap.jl")
include("peak_detect.jl")

NG, M::Int32 = 1, 1000
Hₛ, Tₚ, Tᵢ, Tₑ, γ, g::Float64 = 3.0, 8.0, 0.2, 100.0, 3.3, 9.81

fⱼ,Sⱼ = spectrum(Hₛ,Tₚ,γ,Tᵢ,Tₑ, M) # Generate JONSWAP spectrum for pair (Hₛ,Tₚ) ∈ [Tᵢ,Tₑ]
N = length(fⱼ)

ωₚ = 2π / Tₚ;   kₚ = ωₚ^2/g;    λₚ = 2π / kₚ
ω = 2π * fⱼ;    k = 1/g * ω.^2 

t = range(-10*Tₚ,10*Tₚ,M)
t = collect(t)   
ηₜ = zeros(Float64,M)
# x = range(-5*λₚ,5*λₚ,M);   ηₓ = zeros(Float64,M)

df = abs(fⱼ[end]-fⱼ[1])/(N-1)
η̂ = zeros(Float64,N)
# η̂ = sqrt.(Sⱼ*df/(2*π^2))
# η̂ₘₓ = findmax(η̂)[1]
# A = η̂/η̂ₘₓ

m₀ = Hₛ^2/16
H₀= sqrt(2*m₀*log(N))

Random.seed!(146)
distr = Rayleigh(1) # Rayleigh distribution with unit scale (i.e. Rayleigh(1))
ϕ = π/2 * rand(Uniform(-1,1),NG)
# n = range(0,NG-1)
# ϕ = π/NG * n
A = rand(distr,N)

for i ∈ 1:NG
    for n ∈ 1:N
        η̂[n] = H₀/m₀ * Sⱼ[n]*df  
        for m ∈ 1:M
            ηₜ[m] = ηₜ[m] + A[n]*η̂[n] * cos(-ω[n]*t[m]+ϕ[i])
            # ηₜ[m] = ηₜ[m] + η̂[n] * cos(-ω[n]*t[m]+ϕ[i])
        end
    end
end

tₘ⁺, tₘ¯, ηₘ⁺, ηₘ¯, η̇ₘ, η̈ₘ = peaks(ηₜ,t)

# maxima = [ηₘ⁺; ηₘ¯]
# t_maxima = [tₘ⁺; tₘ¯]

# for i ∈ 1:NG
#     for j ∈ 1:N
#         for m ∈ 1:M
#             η[m] = η[m] + η̂[j] * cos(k[j]*x[m]+ϕ[i])
#         end
#     end
# end

freq, mag, _ = one_side_asp(ηₜ,t)

display(plot(fⱼ, η̂, xlabel = L"f~[Hz]", ylabel = L"\tilde{\eta}~[]", label = L"Wave~Group", linewidth=3))
display(plot(t, ηₜ, xlabel = L"t~[s]", ylabel = L"\eta~[]", label = L"Wave~Group", linewidth=3))
# display(plot(x, ηₓ, xlabel = L"x~[m]", ylabel = L"\eta~[]", label = L"Wave~Group", linewidth=3))
plot(freq, mag, xlabel = L"f~[Hz]", ylabel = L"Magnitude", label = "Single-sided Spectrum", linewidth=3)

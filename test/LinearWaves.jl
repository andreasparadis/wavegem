# Monochromatic waves at certain spatil point (x=0)
using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles, FFTW
using Statistics, Random, Distributions
using Dates

include("wave_theory.jl")
include("signal_processing.jl")

d, H, T::Float64 = 100, 1, 8.0
ρ, g::Float64 = 1025.0, 9.81
κ, υ ::Float64 = 0.0, 0.0
regime ::Int8 = 0

ω = 2π / T

# Determine regime (Deep/Intermediate/Shallow water)
λ₀ = g/(2*π) * T^2 # Deep water wavelength / Initial guess
if d > λ₀/2
    print("Deep water \n")
    λ = λ₀
elseif d>λ₀/20 && d<=λ₀/2 
    print("Intermediate water \n")
    λ,_ = dispersion(T,d)
    regime = 1
else
    print("Shallow water \n")
    λ,_ = dispersion(T,d)
    regime = 2
end

κ = 2π/λ;   υ = ω/κ

# Time discretization 
tₑ = 30*T
# tₑ = round(2*Ldom/υ /10)*10
dt = 0.2
Nₜ = Int64(round(tₑ/dt))+1
t = zeros(Float64,Nₜ)
[t[i] = (i-1)*dt for i ∈ 1:Nₜ]
# Depth discretization
dz = 1
Nz = Int64(round(d/dz))+1
z = zeros(Float64,Nz)
[z[i] = -(i-1)*dz for i ∈ 1:Nz]

θ = ω*t;    csθ = cos.(θ);   snθ = sin.(θ)

# Horizontal spatial discretization
Ldom = 10*round(λ) # [m] Domain length
Cr = 0.5
dx = dt*υ/Cr
Nₓ = Int64(round(Ldom/dx))+1

# Amplitudes
η̂ = H/2 * cosh.(κ*(z.+d))/cosh(κ*d)
Φ̂ = g/ω*H/2 * cosh.(κ*(z.+d))/cosh(κ*d)
û = ω*H/2 * cosh.(κ*(z.+d))/sinh(κ*d)
ŵ = ω*H/2 * sinh.(κ*(z.+d))/sinh(κ*d)
p̂ = ρ*g*H/2 * cosh.(κ*(z.+d))/cosh(κ*d)
Û = mean(û)

# Compute variables
zᵢ = -100 # Pick value z ∈ [0:-d]
zid = Int64(round(-zᵢ/dz)) + 1

ηᵢ = η̂[zid].*snθ # Wave profile at depth z
Φ = Φ̂[zid].*csθ  # Velocity potential
u = û[zid].*snθ  # Horizontal velocity
w = ŵ[zid].*csθ  # Vertical velocity
p = p̂[zid].*snθ  # Dynamic pressure
U = Û .* snθ     # Uniform velocity

η = η̂[1].*snθ    # Surface elevation
Φₜ = -g*η

ηₜ = zeros(Float64, Nₜ)
for i ∈ 2:Nₜ-1
    ηₜ[i] = (η[i+1] - η[i-1]) / (t[i+1]-t[i-1])
end

## PLOTS
plt0 = plot(t, η, xlab = L"t~[s]", ylab = L"\eta~[m]", lab = L"η(z=0, t)", lw=1) 
plot!(t, ηᵢ, lab = L"η(z="*string(zᵢ)*L", t)", lw=1) 
display(plt0)

plt1 = plot(t, Φ, xlab = L"t~[s]", ylab = " ", lab = L"\Phi(z="*string(zᵢ)*L", t)", lw=1) 
plot!(t, Φₜ, lab = L"\frac{\partial \Phi}{\partial t}(z="*string(zᵢ)*L", t)", lw=1)
display(plt1)

plt2 = plot(t, u, xlab = L"t~[s]", ylab = L"Velocity~[m/s]", lab = L"u(z="*string(zᵢ)*L", t)", lw=1) 
plot!(t, w, lab = L"w(z="*string(zᵢ)*L", t)", lw=1) 
display(plt2)

plt3 = plot(t, p, xlab = L"t~[s]", ylab = L"Pressure~[Pa]", lab = L"p_D(z="*string(zᵢ)*L", t)", lw=1) 
display(plt3)
# plt3 = plot(t, p.+1e6.-ρ*g*zᵢ, xlab = L"t~[s]", ylab = L"\Pressure~[Pa]", lab = L"p_D(z="*string(zᵢ)*L", t)", lw=1) 
# display(plt3)

plt4 = plot(η̂, z, xlab = L"Amplitude", ylab = L"z~[m]", lab = L"\hat{\eta}(z)", lw=1)
plot!(û, z, lab = L"\hat{u}(z)", lw=1)
plot!(ŵ, z, lab = L"\hat{w}(z)", lw=1, line=:dash)
display(plt4)

plt5 = plot(p̂/1e6, z, xlab = L"Dynamic~pressure~[\times 10^6 ~Pa]", ylab = L"z~[m]", lab = L"\hat{p}(z)", lw=1)
display(plt5)

plt6 = plot(Φ̂, z, xlab = L"Velocity~Potential~[m^2/s]", ylab = L"z~[m]", lab = L"\hat{\Phi}(z)", lw=1)
display(plt6)

plt7 = plot(t, ηₜ, xlab = L"t~[s]", ylab = L"Velocity~[m/s]", lab = L"\eta_t(z="*string(zᵢ)*L", t)", lw=1) 
display(plt7)

plt8 = plot(t, u*d, xlab = L"t~[s]", ylab = L"Flux~[m^2/s]", lab = L"d\times u(z="*string(zᵢ)*L", t)", lw=1) 
display(plt8)

plt9 = plot(t, U, xlab = L"t~[s]", ylab = L"U~[m/s]", lab = L"U(t)", lw=1) 
display(plt9)

## WRITE OUTPUTS
# Write surface elevation in text file
fid::String = "lin_eta.txt" # File name
open("library/linear/"*fid, "w")
head = ["Monochromatic wave" " "]
info1 = [string(dt)*" "*string(Nₜ)*" 1" " "]
info2 = [string(0.5) " "]
writedlm("library/linear/"*fid, [head; info1; info2; t η], '\t')

# Write horizontal velocity in text file
fid::String = "lin_u.txt" # File name
open("library/linear/"*fid, "w")
writedlm("library/linear/"*fid, [head; info1; info2; t u], '\t')

# Write vertical velocity in text file
fid::String = "lin_w.txt" # File name
open("library/linear/"*fid, "w")
writedlm("library/linear/"*fid, [head; info1; info2; t w], '\t')

# Write potential flux in text file
fid::String = "lin_Phit.txt" # File name
open("library/linear/"*fid, "w")
writedlm("library/linear/"*fid, [head; info1; info2; t Φₜ], '\t')

# Write flux in text file
fid::String = "lin_flux.txt" # File name
open("library/linear/"*fid, "w")
writedlm("library/linear/"*fid, [head; info1; info2; t u*d], '\t')
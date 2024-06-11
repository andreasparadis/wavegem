# Plot the final surface elevation from the OW3D output text file "eta0_irregular"
# Spectrum and envelope of the 2nd half of the signal (based on return period of highest frequency)
using Plots
using DataFrames
using LaTeXStrings
using LinearAlgebra
using DelimitedFiles
using FFTW, FourierAnalysis, DSP
using Interpolations

include("signal_processing.jl")
include("jonswap.jl")

fid::String = "run_8/eta0_irregular.txt"
cont = Array{Float64}(undef,0)
# Hₛ, Tₚ, Tᵢ, Tₑ, γ::Float64 = 3.0, 8.0, 0.4, 1000, 3.3 # Real scale

# fⱼ,Sⱼ = spectrum(Hₛ,Tₚ,γ,Tᵢ,Tₑ) # Generate JONSWAP spectrum for pair (Hₛ,Tₚ) ∈ [Tᵢ,Tₑ]

open("library/"*fid, "r")
cont = readdlm("library/"*fid, '\t', Float64, skipstart=1, '\n')

N = length(cont[:,1])
t = zeros(Float64,N,1)
η = zeros(Float64,N,1)

t = cont[:,1]
η = cont[:,2]
dt = (t[end]-t[1])/(N-1)

if mod(N,2)==1
    print("The signals do not have even number of elements.\n")
    print("For combatibility with Fourier transform the relevant vectors are appended linearly.\n")
    print("The scheme used assumes dt is constant.\n")
    print("-----------------------------------------------------------------------------------.\n")

    push!(t,t[end]+dt)
    push!(η,2*η[end]-η[end-1])
    N = N+1
end

# freq, mag = one_side_asp!(η₀,t)

# Get the half end of the signal (start at one return period of highest frequency) ! REFACTOR !
M = Int(floor(95*N/100))
# To just type x[M:end] in the FFT call below
if mod(M,2)==0
    M = M+1
end

freq, mag,df = one_side_asp!(η[M:end],t[M:end])
ηᴴ = hilbert(η[M:end])
ηᴴₑₙᵥ = abs.(ηᴴ)

m = length(ηᴴₑₙᵥ)
τ = t[M:end]
dη = zeros(Float64,m,1)
ddη = zeros(Float64,m,1)
imax = zeros(Float64,1)
imin = zeros(Float64,1)
τmax = zeros(Float64,1)
ηmax = zeros(Float64,1)
τmin = zeros(Float64,1)
ηmin = zeros(Float64,1)

for i ∈ 2:m-1
    n = M-1+i
    dη[i] = (η[n+1] - η[n-1]) / (τ[i+1]-τ[i-1])
end

dη[1] = 2*dη[2]-dη[3]

for i ∈ 2:m-1
    n = M-1+i
    ddη[i] = (dη[i+1] - dη[i-1]) / (τ[i+1]-τ[i-1])
    if sign(dη[i]) != sign(dη[i-1])
        if ddη[i-1] < 0
            push!(imax,i-1)
            push!(τmax,τ[i-1])
            push!(ηmax,η[n-1])
        else
            push!(imin,i-1)
            push!(τmin,τ[i-1])
            push!(ηmin,η[n-1])
        end
    end
end

ddη[end] = 2*ddη[end-1]-ddη[end-2]

imax = imax[2:end]
imin = imin[2:end]
τmax = τmax[2:end]
ηmax = ηmax[2:end]
τmin = τmin[2:end]
ηmin = ηmin[2:end]


max_itp = interpolate(ηmax, BSpline(Cubic(Free(OnGrid()))))
min_itp = interpolate(ηmin, BSpline(Cubic(Free(OnGrid()))))
# max_eval = max_itp(imax)
# min_eval = min_itp(imin)

res = (max_itp(1:0.1:28)+min_itp(1:0.1:28))/2

display(plot(t[M:end], [η[M:end] ηᴴₑₙᵥ], xlab = L"t~[s]", ylab = L"\eta~[m]", lab = [L"\eta (t)" L"Envelope"]))
display(plot(freq[1:floor(Integer,0.5/df)], mag[1:floor(Integer,0.5/df)], xlab = L"f~[Hz]", ylab = L"Magnitude", lab = "Single-sided Spectrum", lw=3))
plot(τ, η[M:end], xlab = L"t~[s]", ylab = L"\eta~[m]", lab = "η(t)")
plot!(τmax, ηmax, seriestype=:scatter, lab = "maxima")
display(plot!(τmin, ηmin, seriestype=:scatter, lab = "minima"))

plot(max_itp)
plot!(min_itp)
plot!(ηmax, seriestype=:scatter, lab = "maxima")
display(plot!(ηmin, seriestype=:scatter, lab = "minima"))

plot(res)

# display(plot(τ,[η[M:end] dη],xlab = L"t~[s]", ylab = "", lab = [L"\eta" L"\frac{d\eta}{dt}"]))
# plot(τ,[dη ddη],xlab = L"t~[s]", ylab = "", lab = [L"\frac{d\eta}{dt}" L"\frac{d^2\eta}{dt^2}"])
# plot!(fⱼ, Sⱼ, xlab = L"f~[Hz]", ylab = L"S(f)~[m^2 s]", lab = L"JONSWAP", lw=3)
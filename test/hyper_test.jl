using Plots
using DataFrames
using LaTeXStrings
using LinearAlgebra
using FFTW

include("signal_processing.jl")

n::Int64 = 2^10
z::Float64 = 10

k = 0.5

d = range(10,1000,n)

# x = zeros(Float64,n)
t = range(-π,π,n)
f₁ = zeros(Float64,n)
f₂ = zeros(Float64,n)


for i=1:n
    # x[i] = k*(d[i] + z)
    if t[i] ≠ 0
        f₁[i] = sin(4π*t[i]) / (4π*t[i])
    else
        f₁[i] = 1
    end
    # f₁[i] = cosh(x[i])/sinh(k*d[i])
    # f₂[i] = exp(k*z)
end

# freq, mag = one_side_psd(f₁,t)
freq, mag, ϕ₁, freq2, mag2, ϕ₂ = full_fft_psd(f₁,t)

display(plot(t./π, f₁, xlab = L"\frac{t}{\pi}", ylab = "A", lab = L"sinc(t)"))
display(plot(freq, mag, xlab = L"f~[Hz]", ylab = L"Magnitude", lab = "Single-sided Spectrum", lw=3))
display(plot(freq2, mag2, xlab = L"f~[Hz]", ylab = L"Magnitude", lab = "Double-sided Spectrum", lw=3))
plot(freq, ϕ₁, xlab = L"f~[Hz]", ylab = L"Phase", lab = "Double-sided Spectrum", lw=3)
# savefig("sinc.svg")
# plot!(x, f₂, xlabel = "x", ylabel = "y", label = L"e^x")
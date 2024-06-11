using Plots
using DataFrames
using LaTeXStrings
using LinearAlgebra
using DelimitedFiles
using FFTW
using Statistics

include("signal_processing.jl")

open("library/IrregAiry_surface.txt", "r")
cont = readdlm("library/IrregAiry_surface.txt", '\t', Float64, '\n')

t = cont[:,1]
η = cont[:,2]

psd_η = one_side_asp!(η,t)
m₀ = var(η)
Hₘ₀ = 4*sqrt(m₀)

display(plot(t, η, xlabel = L"t~[s]", ylabel = L"\eta~[m]", label = L"Irregular~Waves", linewidth=3))
display(plot(t[1:1000], η[1:1000], xlabel = L"t~[s]", ylabel = L"\eta~[m]", label = L"Irregular Waves", linewidth=3))
plot(psd_η[1], psd_η[2], xlabel = L"f~[Hz]", ylabel = L"Magnitude", label = "Single-sided Spectrum", linewidth=3)
# plot(f2, P2, xlabel = L"f~[Hz]", ylabel = L"Magnitude", label = "Double-sied Spectrum", linewidth=3)
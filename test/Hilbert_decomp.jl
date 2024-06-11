using Plots, LaTeXStrings, DelimitedFiles, Dates, LinearAlgebra
using FFTW, DSP, Statistics, BSplineKit

# Include necessary scripts for functions
include("signal_processing.jl")
include("jonswap.jl")
include("peak_detect.jl")
include("directories.jl")
include("wave_theory.jl")

#############################################################################################
# Import time-domain data for each pi/2 phase shifted realization
𝚶⃗ = Array{Float64}(undef,0,2)
# 1st column: Time    2nd column: Data
A00, A05, A10, A15 = (𝚶⃗[:] for _ = 1:4)
# A00
fid ::String = Decpath*"eta0" # Path to file
open(fid, "r")
A00 = readdlm(fid, '\t', Float64, '\n')
# A05
fid ::String = Decpath*"eta1" # Path to file
open(fid, "r")
A05 = readdlm(fid, '\t', Float64, '\n')
# A10
fid ::String = Decpath*"eta2" # Path to file
open(fid, "r")
A10 = readdlm(fid, '\t', Float64, '\n')
# A15
fid ::String = Decpath*"eta3" # Path to file
open(fid, "r")
A15 = readdlm(fid, '\t', Float64, '\n')

#############################################################################################
ηᴴ = hilbert(A00[:,2]);   ηᴴₑₙᵥ = abs.(ηᴴ)
A05 = A00[:,:];   A10 = A00[:,:];  A15 = A00[:,:] 
A05[:,2] = real(ηᴴₑₙᵥ.*exp.(1im*(angle.(ηᴴ).+π/2)))
A10[:,2] = real(ηᴴₑₙᵥ.*exp.(1im*(angle.(ηᴴ).+π)))
A15[:,2] = real(ηᴴₑₙᵥ.*exp.(1im*(angle.(ηᴴ).+3π/2)))
A05[:,2] .-= mean(A05);    A10[:,2] .-= mean(A10);    A15[:,2] .-= mean(A15)
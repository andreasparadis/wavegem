using Plots, LaTeXStrings, DelimitedFiles, Dates, LinearAlgebra
using FFTW, DSP, Statistics, BSplineKit

gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))

# Include necessary scripts for functions
include("func/text_process.jl")
include("func/signal_processing.jl")
include("func/jonswap.jl")
include("func/wave_theory.jl")
#############################################################################################
# Methods Definition
mutable struct Wave 
    f::Float64
    T::Float64
    ω::Float64
    λ::Float64
    κ::Float64
    υᶜ::Float64
end

# Global Variables
const ρ, g::Float64 = 1025.0, 9.81  # Water density [kg/m³], gravity [m/s²]
const Ldom, d::Float64 = 20.0, 1.0  # [m] Domain length, Water depth

Tₑ, γ, fcut::Float64 =  128, 3.3, 4.0  # JONSWAP parameters
prb_id::Int8 = 4        # Probe id at phase focusing location (Channel No)
prb_pos::Float64 = 8.0  # Probe position [m]

cfld = "HS006TP1TR128"                      # Case folder
# cfld = joinpath("HS006TP1TR128","CORR")
sigf = ("A00.txt", "A05.txt", "A10.txt", "A15.txt") # Signals to decompose

# Paths    
libpath = joinpath(pwd(),"library") # Project's library folder (don't edit)
Parent = joinpath(libpath,"UCL","Waves",cfld)
Decpath = joinpath(Parent,"Decomposition")
CBpath = joinpath(Parent,"CLB")
DecFigs = joinpath(Decpath,"Figures")

#############################################################################################
# Read 1st order decomposition file
fid = joinpath(Decpath,"eta_lin") # Time series
fcont = parse_fxw(fid, 1)
time = fcont[:,1] 
η¹ = fcont[:,2]  
nₜ = length(time)  
tₑ = time[end]-time[1]    # [s] Simulation duration
dt = tₑ/(nₜ-1)

fid = joinpath(Decpath,"eta_lin_spec") # Spectral analysis
fcont = parse_fxw(fid, 1)
FREQ1 = fcont[:,1] 
MAG1 = fcont[:,2]
ARG1 = fcont[:,3]

# itp_sp = interpolate(FREQ1, MAG1, BSplineOrder(4))
itp_phi = interpolate(FREQ1, unwrap(ARG1), BSplineOrder(4))

# JONSWAP spectra
Tᵢ = 1/fcut
long_wave = wave_qnts(Tₑ,d)
ω⁻ = long_wave.ω
ω⁺ = 2π*fcut
dω = 2π/((nₜ-1)*dt)
Nₛ = Int64(round((ω⁺-ω⁻)/dω))

# fₛ = 1/dt  # [Hz] Sampling frequency
# Nₛ = (fₛ*tₑ+1)/4
# fₛ = (4*Nₛ-1)/tₑ  # [Hz] Sampling frequency
# df = 1/tₑ       # Frequency resolution
# fcut = 2^9/Tₑ   # [Hz] Cut-off frequency 

# Target
Hₛᵗ = 0.06
Tₚᵗ = 1.00
fₚᵗ = 1/Tₚᵗ

println("Target spectrum (What output should be like)")
fⱼ,Sⱼᵗ = spectrum(Hₛᵗ,Tₚᵗ,γ,Tᵢ,Tₑ,Nₛ) # JONSWAP target spectrum
dfⱼ = abs(fⱼ[end]-fⱼ[1])/(length(fⱼ)-1)
η̂ᵗ = sqrt.(2*Sⱼᵗ*dfⱼ)

# Out
# Hₛᵒ = 4*std(η¹)
# fₚᵒ = FREQ1[findmax(MAG1)[2]]
# Tₚᵒ = 1/fₚᵒ

fid = joinpath(Decpath,"JONSWAP_pars") # Output JONSWAP parameters
fcont = parse_fxw(fid, 1)
Hₛᵒ, Tₚᵒ = fcont[1], fcont[2]
fₚᵒ = 1/Tₚᵒ
ωₚᵒ = 2π*fₚᵒ

println("Output Spectrum (1st order elevation - latest run)")
_,Sⱼᵒ = spectrum(Hₛᵒ,Tₚᵒ,γ,Tᵢ,Tₑ,Nₛ) # output JONSWAP interpolation
η̂ᵒ = sqrt.(2*Sⱼᵒ*dfⱼ)

# Amplitude Correction
af_tgt = η̂ᵗ[:]
af_in = af_tgt[:]
af_out = η̂ᵒ[:]
# af_out = itp_sp.(fⱼ)

af_in = af_in .* af_tgt ./ af_out

for i ∈ 1:Nₛ 
    if isequal(af_in[i], NaN)
        af_in[i] = 0
    end
end

fₚⁱ = fⱼ[findmax(af_in)[2]]
ωₚⁱ = 2π*fₚⁱ
Sⱼⁱmax = maximum(af_in) # =ag^2/ωₚ^5 * exp(-5/4) * γ

α = Sⱼⁱmax * ωₚⁱ^5 / (g^2 * exp(-5/4)*γ)
Hₛⁱ = 2.227* sqrt(α) * g / ωₚⁱ^2

# Phase Correction (for Focused WGs)
phi_tgt = 0*zeros(Float64,Nₛ) # 0, π/2, π, 3π/2
phi_in = phi_tgt[:]
phi_out = itp_phi.(fⱼ)

phi_in = phi_in .- (phi_tgt .- phi_out)

#  Calculate phase difference between focus location and X = 0
phi_diff = -2π*prb_pos/g * fⱼ

#  Calculate phases at X = 0 by subtracting phase difference
phi_in = phi_in .- phi_diff

println("Input spectrum (What should be passed to the wavemaker)")
_,Sint = spectrum(0.069,1/fₚⁱ,γ,Tᵢ,Tₑ,Nₛ)
η̂int = sqrt.(2*Sint*dfⱼ)

plt_sp = plot(fⱼ,af_tgt, lab="Target", lw=2)
plot!(fⱼ,af_in, lab="Input", lw=2)
plot!(fⱼ,η̂int, lab="Fit Input", ls=:dash, lw=1)
plot!(fⱼ,af_out, lab="Output", lw=2)

plt_arg = plot(fⱼ,phi_tgt, lab="Target", lw=2)
plot!(fⱼ,phi_in, lab="Input", lw=2)
plot!(fⱼ,phi_out, lab="Output", lw=2)

display(plt_sp)
display(plt_arg)
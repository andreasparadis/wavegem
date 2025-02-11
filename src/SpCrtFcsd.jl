using Plots, LaTeXStrings, DelimitedFiles, Dates, LinearAlgebra
using FFTW, DSP, Statistics, BSplineKit
using CurveFit: curve_fit, linear_fit, Polynomial, expsum_fit, ExpFit

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

cfld = joinpath("HS006TP1TR128","Events","EV1","FCSD")      # Case folder
sigf = ("A00.txt", "A05.txt", "A10.txt", "A15.txt") # Signals to decompose

# Paths    
libpath = joinpath(pwd(),"library") # Project's library folder (don't edit)
Parent = joinpath(libpath,"UCL","Waves",cfld)
Decpath = joinpath(Parent,"Decomposition")
CBpath = joinpath(Parent,"CLB")
DecFigs = joinpath(Decpath,"Figures")
FRpath = joinpath(libpath,"UCL","Waves","HS006TP1TR128","Random","Decomposition","ReEvent")
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

itp_sp = interpolate(FREQ1, MAG1, BSplineOrder(4))
itp_phi = interpolate(FREQ1, unwrap(ARG1), BSplineOrder(4))

# Spectral parameters
Tᵢ = 1/fcut
long_wave = wave_qnts(Tₑ,d)
ω⁻ = long_wave.ω
ω⁺ = 2π*fcut
dω = 2π/((nₜ-1)*dt)
Nₛ = Int64(round((ω⁺-ω⁻)/dω))+1
fⱼ = range(1/Tₑ,1/Tᵢ,Nₛ)

# fₛ = 1/dt  # [Hz] Sampling frequency
# Nₛ = Int64(round((fₛ*tₑ+1)/4))
# fₛ = (4*Nₛ-1)/tₑ  # [Hz] Sampling frequency
# df = 1/tₑ       # Frequency resolution
# fcut = 2^9/Tₑ   # [Hz] Cut-off frequency 
# fⱼ = range(1/Tₑ,1/Tᵢ,Nₛ)

# Target
println("Target spectrum (What output should be like)")
fid = joinpath(FRpath,"EV1_FCSD_FR.fronts")
fcont = parse_fxw(fid, 3)
FREQ = fcont[:,1] 
MAG = fcont[:,2]
ARG = fcont[:,4]
itp_spt = interpolate(FREQ, MAG, BSplineOrder(4))
itp_phit = interpolate(FREQ, unwrap(ARG), BSplineOrder(4))

# Out
Hₛᵒ = 4*std(η¹)
fₚᵒ = FREQ1[findmax(MAG1)[2]]
Tₚᵒ = 1/fₚᵒ

# Amplitude Correction
af_tgt = itp_spt.(fⱼ)
af_in = af_tgt[:]
af_out = itp_sp.(fⱼ)
af_out = low_pass_filter(af_out,4,100,2^8)

af_in = af_in .* af_tgt ./ af_out

for i ∈ 1:Nₛ 
    if isequal(af_in[i], NaN)
        af_in[i] = 0
    end
end

fₚⁱ = fⱼ[findmax(af_in)[2]]
ωₚⁱ = 2π*fₚⁱ
Sⱼⁱmax = maximum(af_in)

# Phase Correction (for Focused WGs)
phi_tgt = 0*zeros(Float64,Nₛ) # 0, π/2, π, 3π/2
phi_in = phi_tgt[:]
phi_out = itp_phi.(fⱼ)

phi_in = phi_in .- (phi_tgt .- phi_out)

#  Calculate phase difference between focus location and X = 0
phi_diff = -2π*prb_pos/g * fⱼ

#  Calculate phases at X = 0 by subtracting phase difference
phi_in = phi_in .- phi_diff

H2 = Nₛ*af_in.*exp.(1im*phi_in)
H2 = vcat(H2,0,conj(H2[end:-1:2]))
ETA = real(ifft(H2))

plt_sp = plot(fⱼ,af_tgt, lab="Target", lw=2)
plot!(fⱼ,af_in, lab="Input", lw=2)
plot!(fⱼ,af_out, lab="Output", lw=2)

plt_arg = plot(fⱼ,phi_tgt, lab="Target", lw=2)
plot!(fⱼ,phi_in, lab="Input", lw=2)
plot!(fⱼ,phi_out, lab="Output", lw=2)

plt_out = plot(xlab=L"t~[s]", ylab=L"\eta~[m]")
plot!(range(0,tₑ,length(ETA)),ETA, lab="In")

display(plt_sp)
display(plt_arg)
display(plt_out)

suffix = "EV1_FCSD_FR_CORR_1"
fid = joinpath(Decpath,"EV1",suffix*".fronts")
open(fid, "w")
head0 = [suffix "" "" ""]
head1 = ["f" "a" "angle" "phase"]
head2 = ["Hz" "m" "rad" "rad"]
wcont = [head0; head1; head2; fⱼ af_in zeros(Float64,Nₛ) phi_in]
writedlm(fid, wcont, '\t')
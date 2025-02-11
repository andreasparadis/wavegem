using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles
using FFTW, DSP, Statistics, BSplineKit
import ColorSchemes.darkrainbow

gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))
cb = darkrainbow

# Include necessary scripts for functions
include("text_process.jl"); include("signal_processing.jl")

# Paths for probe measurements
fdate = "2024-10-09"
dname = "EV1_FCSD_CORR_00"

parent = joinpath(pwd(),"library","UCL","Model")
prbpath = joinpath(parent,"probes",fdate,dname*".txt")
prbfigs = joinpath(parent,"probes",fdate)

###########################################################################
# Read measurements from wave probes
fₚ = 1 # [Hz] Wave peak frequency

prb_cont = parse_fxw(prbpath, 1)

tₚᵣ = prb_cont[:,1]
tₚᵣ = tₚᵣ .- tₚᵣ[1]

Vres = 0.02 # [m/Volt]
Zₚᵣ = Vres * prb_cont[:,6]

###########################################################################
# Spectral analysis
FR, MAG, phi, df, T, Nfft, H = one_side_asp(Zₚᵣ, tₚᵣ)
MAGₚ = findmax(MAG)[1]

# Alternative processing / Re-sampling
tₑ = 128.0    # [s] Simulation duration
dtᵢ = tₚᵣ[2] - tₚᵣ[1]
fₛᵢ = 1/dtᵢ      # [Hz] Original sampling frequency
Nf = nextpow(2,Int64(round(tₑ*fₛᵢ+1)))    # No of frequency components

df = 1/tₑ       # Frequency resolution
fₛ = (Nf-1)/tₑ  # [Hz] Sampling frequency
fᶜ = 5          # [Hz] Cut-off frequency

dt = 1/fₛ
L = Int64(round(tₑ/dt))+1   # Length of the time signal
t = zeros(Float64,L)
[t[i] = (i-1)*dt for i ∈ 1:L]

## Low-pass filter
# fₛˡᵖ = 4*fᶜ # Filter sampling frequency
fₛˡᵖ = fₛ/2 # Filter sampling frequency
Ntaps = 2^5+1   # No of taps
Zf = low_pass_filter(Zₚᵣ,fᶜ,fₛˡᵖ,Ntaps)

## B-spline interpolation of signal 
itp_xᵢ = interpolate(tₚᵣ, Zf, BSplineOrder(4))
Z = itp_xᵢ.(t)

FRI, MAGI, phiI, _, _, NfftI, _ = one_side_asp(Z, t)
###########################################################################
# Plots
## Probe 4
plt_prb4 = plot(xlab="t [s]", ylab="η [m]", title = "Surface elevation", palette=[cb[11];cb[4];cb[8]])
plot!(tₚᵣ,Zₚᵣ,lw=2,lab="Probe #4")
plot!(t,Z,lw=2,lab="Re-sampled")
plot!(xlim=(40,48))
display(plt_prb4)

# spectrum
plt_SPZ = plot(xlab=L"f~[Hz]", ylab=L"Magnitude~[m]", palette=[cb[11];cb[4];cb[8]])
plot!(FR,MAG,label="Surface elevation", lw=2)
plot!(FRI,MAGI,label="Re-sampled", lw=2)
plot!([fₚ;fₚ],[0;MAGₚ],lab=L"f_p",line=:dashdot, lw=2)
plot!(xlim=(0,10))
display(plt_SPZ)

plt_PHZ = plot(xlab=L"f~[Hz]", ylab=L"Phase~[rad]", palette=[cb[11];cb[4];cb[8]])
plot!(FR,unwrap(phi),label="Original", lw=2)
plot!(FRI,unwrap(phiI),label="Re-sampled", lw=2)
plot!([fₚ;fₚ],[0;MAGₚ],lab=L"f_p",line=:dashdot, lw=2)
plot!(xlim=(0,10))
using Plots, LaTeXStrings, DelimitedFiles, Dates, LinearAlgebra
using FFTW, DSP, Statistics, BSplineKit

gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))

# Include necessary scripts for functions
include("func/text_process.jl")
include("func/signal_processing.jl")
include("func/jonswap.jl")
include("func/peak_detect.jl")
include("func/directories.jl")
include("func/wave_theory.jl")
include("func/event_proc.jl")

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
## Significant wave height [m], Peak period [s], peakedness [-], Cut-off frequency [Hz]
Hₛ, Tₚ, γ, fcut::Float64 = 0.081, 1, 3.3, 4  # JONSWAP parameters 
const ρ, g::Float64 = 1025.0, 9.81  # Sea water density [kg/m³], Accel. of gravity [m/s²]
Ldom, d::Float64 = 10.0, 0.4    # [m] Domain length, Water depth


Decpath =  "/home/andreasp/WAVEGEM/library/JFM/FULL/"
DecFigs = Decpath

tstart = 0; tend = 512
frec, fplot = Bool(1), Bool(1)
sigf = "1.txt"

# Reference signal
fid = Decpath*sigf # Path to file
A00 = parse_fxw(fid, 0)

# Signal truncation
## Truncate t vector
dt = (A00[end,1]-A00[1,1])/(length(A00[:,1])+1)    # t step
ibeg = Int64(round(tstart/dt))+1        # lower limit
iend = Int64(round(tend/dt))            # upper limit
t00 = A00[ibeg:iend, 1] .- A00[ibeg, 1] # Set first value as t zero
tOG = round.(t00*100)/100               # Global t vector for all components
nₜ = length(tOG)

## Truncate vectors of values
A̅00 = A00[ibeg:iend, 2]

# Spectral analysis via FFT
FR1, Â₀₀, ϕ₀₀, _, Tₛ, L = one_side_asp(A̅00, tOG)
df = (FR1[end]-FR1[1])/(L/2)

#############################################################################################
# Isolate, shift around 0 and plot 'pulses' based on peaks
## Find peaks and sort them in descending order
## Isolate peak sections given a specified t range
t_range = 10  # [s]
iNP = 10
shift = 2
res_f = 1

MinPeakVal = 3*std(A̅00)
MinPeakDist = 1/fcut

tₑᵥⁱ, ev_int, tₑᵥ, events, rng, Tₛᵉᵛ, PeakId, PeakPos, PeakVal, NP = 
    ev_identify(A̅00, tOG, MinPeakVal, MinPeakDist, t_range, iNP, shift, res_f)


## Calculate mean pulse
tmean = tₑᵥ[:,NP]
NLN_mean = mean!(ones(2*rng+1,1), events)

#############################################################################################
# Compare signal section and calculate statistical measures
sslen = t_range/2  # Signal section length[s]

ρᴾ, R̅, σᴿ = sect_corr(sslen, dt, events)

#############################################################################################
# Examine most extreme event
Aₑᵥ = ev_int[:,1]
tev = tₑᵥⁱ[:,1]


freq, mag, _,_,_,Nfft, H = one_side_asp(Aₑᵥ, tev)
# flb = 1; fub = Int(round(fcut/df))
# freq = freq[flb:fub]
# mag = mag[flb:fub]

itp_sp = interpolate(freq, mag, BSplineOrder(4))
FR = range(0,fcut, Nfft)
MAG = itp_sp.(FR)
EV = ifft(H)

fid_ev::String = "/home/andreasp/WAVEGEM/library/JFM/EWG/max_event" # File name
open(fid_ev, "w")
writedlm(fid_ev, [tev Aₑᵥ], '\t')

# 4: Peaks on top of original signal
plt_peaks = plot(tOG, A̅00, label = "Non-linear")
plot!(PeakPos, PeakVal, seriestype=:scatter, ms=2, mc=:red, lab = "Peaks",legend=:bottomleft, legendcolumns=2)
plot!(title = "Peaks of Non-linear signal", xlab = "t [sec]", ylab = "η [m]")

# 6: Peak sections and mean of sections
plta =  plot(tₑᵥ, events, legend = false, grid=true)
plot!(tmean,NLN_mean,line = :dash, color = :red, linewidth = 2, label = "Mean")
plot!(title = "Non-linear signal", xlab = "t [sec]", ylab = "η [m]")


plt_spitp = plot(freq, mag, lw=2, xlab = "f [Hz]", ylab = L"Amplitude", lab="Max event")
plot!(FR, MAG, line =:dashdot, lw=2, lab="Interpolation")
plot!(xlim=(0,fcut))


# 3: Interpolated event against original event
lb = PeakId[1]-rng; ub = lb+2*rng 
tpart = tOG[lb:ub]
plt_itp = plot(tev, Aₑᵥ, line = :solid, lw = 2, lab = "Interpolation", legend=:bottomleft)
plot!(twiny(), tpart, A̅00[lb:ub], line =:dash, color =:red, lw = 2, ylab = "η [m]",  lab = "Non-linear", legend=:bottomright)
plot!(title = "Event at Highest Peak: Original vs Interpolated", xlab = "t [sec]")

# 9: 
h = PeakVal[1]
plt_groups = plot(tOG, A̅00, xlab = "t [sec]", ylab = "η [m]",  label = "Non-linear")
for i ∈ 1:min(20, length(PeakVal))
    c = PeakPos[i]
    w = t_range/2
    plot!(Shape([c-w,c-w,c+w,c+w],[-h,h,h,-h]), color=:red, opacity=0.2, label = "WG$i")
end
plot!(legend = false)

# 11:
c = PeakPos[1]; w = t_range
plt_maxev = plot(tOG, A̅00, xlab = "t [sec]", ylab = "η [m]",  label = "Non-linear")
plot!(Shape([c-w,c-w,c+w,c+w],[-h,h,h,-h]), color=:red, opacity=0.3, label = "Event")
# plot!(xlim=(c-10*w, c+10*w))

display(plt_peaks)
display(plta)
display(plt_spitp)
display(plt_itp)
display(plt_groups)
display(plt_maxev)
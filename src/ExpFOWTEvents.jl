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

γ, fcut::Float64 = 3.3, 4  # JONSWAP parameters 
mrk_id::Int8 = 1        # Marker ID
prb_pos::Float64 = 8.0  # Probe position [m]
tdur::Float64 = 1200.0  # Signal time range to consider [s]
evID::Int8 = 5          # id of event to save (desc. order from heighest peak)

# Data locations - Edit accordingly
cfld::String = joinpath("Model","cam","2024-10-10")   # Case folder (Focused Waves)
sigf = ("HS006TP1TR128_00-X.txt", "HS006TP1TR128_05-X.txt", 
        "HS006TP1TR128_10-X.txt", "HS006TP1TR128_15-X.txt") # Signals to decompose

libpath = joinpath(pwd(),"library","UCL") # Library folder
Parent = joinpath(libpath,cfld)
evdir = joinpath(Parent,"Surge_EV$evID")

if !isdir(evdir)
    mkdir(evdir)
end

# Flags
tshift = true  # True for experimental data

#############################################################################################
# Import t-domain data for each pi/2 phase shifted realization
## Reference signal
fid = joinpath(Parent,sigf[1]) # Path to file
fcont = parse_fxw(fid, 0)
t00 = fcont[:,1]
A00 = fcont[:,mrk_id+1]

## Phase shifted signals
∅ = zeros(Float64, length(A00[:,1]),3);     A = ∅[:,:]
for i ∈ 1:3
    fid2 = joinpath(Parent,sigf[i+1])
    cont = parse_fxw(fid2, 0)
    A[:,i] = cont[:,mrk_id+1] 
end
A05 = A[:,1];   A10 = A[:,2];   A15 = A[:,3]

##########################################################################################
# Signal truncation
## Truncate t vector
Tᵢ = 1/fcut             # Cut-off period of spectrum [s]
λ⁻ = g/(2π) * Tᵢ^2     # Shortest wavelength (Deep water) [m]
υ⁻ = λ⁻/Tᵢ              # Slowest wave [m/s]
tstart = prb_pos/υ⁻; tend = tstart+tdur

dt = t00[2]-t00[1]                      # t step
ibeg = Int64(round(tstart/dt))+1        # lower limit
iend = Int64(round(tend/dt))            # upper limit
tOG = t00[ibeg:iend] .- t00[ibeg]       # Set first value as t zero
nₜ = length(tOG)

## Truncate vectors of values and make them relative
A̅00 = A00[ibeg:iend] .-A00[1]
A̅05 = A05[ibeg:iend] .-A05[1]
A̅10 = A10[ibeg:iend] .-A10[1]
A̅15 = A15[ibeg:iend] .-A15[1]

#############################################################################################
# 3rd order decomposition of original signal (signal 0)
SPEC0, AMPS, ARG, sig_comps, sig_recon, L, fₚ, pltϕ = decomp_3rd(tOG, A̅00.-mean(A̅00), A̅05, A̅10, A̅15)

FR1, Â₀₀, ϕ₀₀ = SPEC0[:,1], SPEC0[:,2], SPEC0[:,3] # Spectral info of signal 0
df = (FR1[end]-FR1[1])/(L/2)

# Event Identification
## Find peaks and sort them in descending order
## Isolate peak sections given a specified t range
t_range = 20  # [s]
iNP = 5
shift = 2
res_f = 1

cstd = 3
MinPeakDist = 1/fcut
MinPeakVal = cstd*std(A̅00)+mean(A̅00)

tₑᵥⁱ, ev_int, tₑᵥ, events, rng, Tₛᵉᵛ, PeakId, PeakPos, PeakVal, NP = 
    ev_identify(A̅00, tOG, MinPeakVal, MinPeakDist, t_range, iNP, shift, res_f)

###################################################################################
# Event limits
lb = PeakId[evID]-rng; ub = lb+2*rng 
tpart = tOG[lb:ub]; Lₚ = length(tpart)

# Examine extreme event
Aₑₑ = A̅00[lb:ub]
tₑₑ = tOG[lb:ub] .- tOG[lb]

# Aₑₑ = ev_int[:,evID]
# tₑₑ = tₑᵥⁱ[:,evID]

fr_ev, mag_ev, phi,_,_,Nfft, H = one_side_asp(Aₑₑ, tₑₑ)
itp_sp = interpolate(fr_ev, mag_ev, BSplineOrder(4))
itp_phi = interpolate(fr_ev, unwrap(phi), BSplineOrder(4))
FR = range(0,fcut, Nfft)
MAG = itp_sp.(FR)
PHI = itp_phi.(FR)
EV = ifft(H)
# dtEV = (tₑₑ[end]-tₑₑ[1])/(Nfft)
tEV = zeros(Float64,Nfft)
tEV[1:Lₚ] = tₑₑ[:]
# [tEV[i] = i*dt for i ∈ Lₚ:Nfft]
icut = Int(round(fcut/(FR[2]-FR[1])))

#############################################################################################
## PLOTS
# 1: Event plot
plt1 = plot(tpart, A̅00[lb:ub], lw = 2, lab = "Non-linear")
plot!(title = "Event$evID", xlab = "t [sec]", ylab = "Surge [mm]")

# 2: Spectral Analysis
flb = 1; fub = Int(round(fcut/df))
plt_sp = plot(FR1[flb:fub], Â₀₀[flb:fub], line =:dashdot, lab="Non-linear")
plot!(xlab = "f [Hz]", ylab =  "Amplitude [mm]")
plot!(xlim=(0,fcut), minorgrid=true)

# 2: Peaks on top of original signal
plt_peaks = plot(tOG, A̅00, label = "Non-linear")
plot!(PeakPos[1:NP], PeakVal[1:NP], seriestype=:scatter, ms=2, mc=:red, lab = "Peaks",legend=:bottomleft, legendcolumns=3)
plot!([tOG[1];tOG[end]],[MinPeakVal;MinPeakVal], lab="$cstd η", mc=:red, ls=:dash)
plot!(title = "Peaks of Non-linear signal", xlab = "t [sec]", ylab = "Surge [mm]")

# 3: Compare event spectrum with full spectrum
fr_tmp, mag_tmp, _ = one_side_asp(Aₑₑ,tₑₑ)
plt_evsp = plot(FR1[flb:fub], Â₀₀[flb:fub], line =:dashdot, lab="Non-linear")
plot!(fr_tmp, mag_tmp, line =:solid, lab="Event")
plot!(xlim=(0,fcut), minorgrid=true)
plot!(title = "Comparison of Spectra", xlab = "f [Hz]", ylab = "Amplitude [mm]")

# 4: Illustrate position of events in the full signal
h = PeakVal[evID]
plt_groups = plot(tOG, A̅00, xlab = "t [sec]", ylab = "Surge [mm]",  label = "Surge")
for i ∈ 1:min(20, length(PeakVal))
    local c = PeakPos[i]
    local w = 10
    plot!(Shape([c-w,c-w,c+w,c+w],[-h,h,h,-h]), color=:red, opacity=0.2)
end
plot!(legend = false)

# 5: Illustrate position of selected event in the full signal
c = PeakPos[evID]; w = 10
plt_maxev = plot(tOG, A̅00, xlab = "t [sec]", ylab = "Surge [mm]",  label = "Surface elevation")
plot!(Shape([c-w,c-w,c+w,c+w],[-h,h,h,-h]), color=:red, opacity=0.3, label = "Event")
plot!(xlim=(max(tOG[1],c-10*w), min(tOG[end],c+10*w)))

# 7: Interpolated event against original event
plt_itp = plot(tEV[1:Lₚ], real(EV[1:Lₚ]), line = :solid, lw = 1, legend=:topleft, lab = "Interpolation")
plot!(twiny(), tpart, A̅00[lb:ub], line =:dash, color =:red, lw = 1, ylab = "Surge [mm]", legend=:topright, lab = "Non-linear")
plot!(title = "Event at Highest Peak: Original vs Interpolated", xlab = "t [sec]")

# 8: Selected event & interpolation - Spectral analysis
plt_mag_itp = plot(fr_ev, mag_ev, lw=2, ylab = L"Amplitude [mm]", lab="Event")
plot!(FR, MAG, line =:dot, lw=2, lab="Interpolation")
plot!(xlim=(0,fcut))

plt_phi_itp = plot(xlab = "f [Hz]", ylab = "ϕ [rad]")
plot!(fr_ev, unwrap(phi), lw=2, lab="Event")
plot!(FR, PHI, lw=2, lab="Interpolation")
plot!(xlim=(0,fcut))
plot!(ylim=(minimum(PHI[1:icut]),maximum(PHI[1:icut])))

plt_spitp = plot(plt_mag_itp, plt_phi_itp, layout = @layout [a; b])

display(plt1)
display(plt_sp)
display(plt_peaks)
display(plt_evsp)
display(plt_groups)
display(plt_maxev)
display(plt_itp)
display(plt_spitp)

# Figures to save
fgnm = joinpath(evdir,"ev_spec")
savefig(plt_sp, fgnm*".png"); savefig(plt_sp, fgnm*".svg")
fgnm = joinpath(evdir,"ev_elev")
savefig(plt1, fgnm*".png");   savefig(plt1, fgnm*".svg")
fgnm = joinpath(evdir,"elev_peaks_nln")
savefig(plt_peaks, fgnm*".png"); savefig(plt_peaks, fgnm*".svg")
fgnm = joinpath(evdir,"events")
savefig(plt_groups, fgnm*".png");   savefig(plt_groups, fgnm*".svg")
fgnm = joinpath(evdir,"event_$evID")
savefig(plt_maxev, fgnm*".png");    savefig(plt_maxev, fgnm*".svg")
fgnm = joinpath(evdir,"sp_ev_$evID")
savefig(plt_spitp, fgnm*".png");    savefig(plt_spitp, fgnm*".svg")

# Write output files
# Where to find this event in the exciting sea-state?
# Only need start and end times.
tᵉₛ = zeros(Float64,NP)
tᵉₑ = zeros(Float64,NP)

for i ∈ 1:NP
    tᵉₛ[i] = tOG[PeakId[i]-rng]
    tᵉₑ[i] = tOG[PeakId[i]+rng]
end

fid = joinpath(evdir,"EV_lim") # File name
open(fid, "w")
head = ["EV" "tstart" "tend"]
writedlm(fid, [collect(range(1,NP)) tᵉₛ tᵉₑ], '\t')
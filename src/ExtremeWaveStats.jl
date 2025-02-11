module ExtremeWaveStats
using Plots, LaTeXStrings, DelimitedFiles, Dates, LinearAlgebra
using FFTW, DSP, Statistics, BSplineKit
import ColorSchemes.darkrainbow

gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))
cb = darkrainbow

include("WAVEGEM.jl")
import .WAVEGEM

# Include necessary scripts for functions
include("func/signal_processing.jl")
include("func/jonswap.jl")
include("func/peak_detect.jl")
include("func/directories.jl")
include("func/wave_theory.jl")
include("func/event_proc.jl")

#############################################################################################
# Import Global Variables & Module specific inputs
ρ, g, Ldom, d, _, _, γ, fcut = WAVEGEM.GlobInp0
pdir, run_id, phi_id, prb_id = WAVEGEM.GlobInp1
case_id, _, rundir = WAVEGEM.GlobPaths[2:4]
Decpath, _, DecFigs = WAVEGEM.GlobPaths[9:11]
OFASTpath = WAVEGEM.GlobPaths[end]
Wave = WAVEGEM.Wave
tstart = WAVEGEM.tstart; tend = WAVEGEM.tend
frec, fplot = WAVEGEM.Dflags 
sigf = WAVEGEM.DecSigs
tshift = false  # True for experimental data

#############################################################################################
# Make directories
evID = 1

evdir = joinpath(Decpath,"EV$evID")
# evdir = joinpath(Decpath,"MaxFair_7")
# evdir = joinpath(Decpath,"MaxCoM")
# evdir = joinpath(Decpath,"MaxPitch")
# evdir = joinpath(Decpath,"MaxWave")

make_dirs(2, Decpath, DecFigs)
make_dirs(2, joinpath(OFASTpath,"ExtElev","$(case_id)"), joinpath(OFASTpath,"ExtElev","$(case_id)","$(run_id)"))

#############################################################################################
# Import t-domain data for each pi/2 phase shifted realization
## Reference signal
fid = joinpath(Decpath,sigf[1]) # Path to file
open(fid, "r")
A00 = readdlm(fid, '\t', Float64, '\n')

## Phase shifted signals
∅ = zeros(Float64, length(A00[:,1]),3);     A = ∅[:,:]
for i ∈ 1:3
    fid = joinpath(Decpath,sigf[i+1]) # Path to file
    open(fid, "r")
    cont = readdlm(fid, '\t', Float64, '\n')
    A[:,i] = cont[:,2]
end
A05 = A[:,1];   A10 = A[:,2];   A15 = A[:,3] 

#############################################################################################
# Signal truncation
## Truncate t vector
dt = A00[2,1]-A00[1,1]                  # t step
ibeg = Int64(round(tstart/dt))+1        # lower limit
iend = Int64(round(tend/dt))            # upper limit
t00 = A00[ibeg:iend, 1] .- A00[ibeg, 1] # Set first value as t zero
tOG = round.(t00*1000)/1000               # Global t vector for all components
nₜ = length(tOG)

## Truncate vectors of values
A̅00 = A00[ibeg:iend, 2]
A̅05 = A05[ibeg:iend]
A̅10 = A10[ibeg:iend]
A̅15 = A15[ibeg:iend]

#############################################################################################
# 3rd order decomposition of original signal (signal 0)
SPEC0, AMPS, ARGS, sig_comps, sig_recon, L, fₚ, pltϕ = decomp_3rd(tOG, A̅00, A̅05, A̅10, A̅15)

FR1, Â₀₀, ϕ₀₀ = SPEC0[:,1], SPEC0[:,2], SPEC0[:,3] # Spectral info of signal 0
df = (FR1[end]-FR1[1])/(L/2)
Val_LNR = sig_comps[:,1]    # 1st order component

T₂₋ = 1/FR1[findmax(AMPS[:,2])[2]] # 2nd- component peak period
T₂₊ = 1/FR1[findmax(AMPS[:,3])[2]] # 2nd+ component peak period

# Instantaneous frequency
𝓗 = hilbert(Val_LNR)     
u = abs.(𝓗)
θ = angle.(𝓗);  θ = unwrap(θ)

omH = zeros(Float64, nₜ)
tom = zeros(Float64, nₜ)
for i ∈ 2:nₜ-1
    omH[i] = (θ[i+1] - θ[i-1]) / (tOG[i+1]-tOG[i-1])
    if omH[i] > 0 
        tom[i] = maximum(A̅00)
    else
        tom[i] = 0
    end
end

ω̅ᵢₙₛₜ = mean(omH)
idᵤₚ = findall(diff(sign.(omH)) .== 2) # Up-crossings of instantaneous frequency
Nᵤₚ = length(idᵤₚ)
ΔTᵉ = diff(tOG[idᵤₚ]);          ΔTᵉ = [tOG[idᵤₚ[1]]; ΔTᵉ]
tΔT = zeros(Float64,Nᵤₚ);       tΔT[1] = tOG[idᵤₚ[1]]/2
evAmax = zeros(Float64,Nᵤₚ);    evAmax[1] = maximum(A̅00[1:idᵤₚ[1]])
for i ∈ 2:Nᵤₚ
    tΔT[i] = (tOG[idᵤₚ[i-1]] + tOG[idᵤₚ[i]])/2
    evAmax[i] = maximum(A̅00[idᵤₚ[i-1]:idᵤₚ[i]])
end

#############################################################################################
# Find peaks above threshold
MinPeakDist = WAVEGEM.Tᵢ
## Non-linear surface eleveation
MinPeakVal = 2*std(A̅00) # Amplitude constraint (Hₛ/2)
PeakVal, PeakPos, PeakId = peaks_max_ext(A̅00, tOG, MinPeakVal, MinPeakDist)
## Linear surface eleveation
MinPeakVal = 2*std(Val_LNR) # Amplitude constraint (Hₛ/2)
PeakVal, PeakPos, PeakId = peaks_max_ext(Val_LNR, tOG, MinPeakVal, MinPeakDist)

#############################################################################################
# Stochastic Processing of Extreme Wave Events
## Isolate, shift around 0 and plot 'pulses' based on peaks
## Sort peaks in descending order and isolate peak sections given a specified t range
t_range = 100  # [s]
iNP = 10
shift = 0   # 0: shift peak at t=0 | 1: start at t=0 | 2: no shift

tₑᵥ, events, rng, PeakId, PeakPos, PeakVal, NP = 
    ev_identify!(A̅00, tOG, PeakVal, PeakPos, PeakId, t_range, iNP, shift)

Ltₑᵥ, Levents, rng, L_PeakId, L_PeakPos, L_PeakVal, LNP = 
    ev_identify!(Val_LNR, tOG, L_PeakVal, L_PeakPos, L_PeakId, t_range, iNP, shift)

## Calculate mean pulse
tmean = tₑᵥ[:,NP]
NLN_mean = mean!(ones(2*rng+1,1), events)

tmeanLNR = Ltₑᵥ[:,LNP]
LNR_mean = mean!(ones(2*rng+1,1), Levents)

#############################################################################################
# Compare signal section and calculate statistical measures
sslen = 5  # Signal section length[s]

ρᴾ, R̅, σᴿ = sect_corr(sslen, dt, Levents)

#############################################################################################
# Examine extreme event
if !isdir(evdir)
    mkdir(evdir)
end

## Event limits
lb = PeakId[evID]
ub = PeakId[evID]

while sign(omH[lb]) == sign(omH[lb-1])
    global lb = lb - 1
end

while sign(omH[ub]) == sign(omH[ub+1])
    global ub = ub + 1
end
# Event exact duration
ΔTₑᵥ = tOG[ub] - tOG[lb]

tpart = tOG[lb:ub]; Lₚ = length(tpart)

Aₑₑ = A̅00[lb:ub]
Aₑₑ¹ = Val_LNR[lb:ub]
tₑₑ = tOG[lb:ub] .- tOG[lb]

fr_ev, mag_ev, phi,_,_,Nfft, H = one_side_asp(Aₑₑ¹, tₑₑ)

#############################################################################################
## PLOTS 
if fplot
    # 1: Peak sections and mean of sections
    plta =  plot(tₑᵥ, events, legend = false, grid=true)
    plot!(tmean,NLN_mean,line = :dash, color = :red, linewidth = 2, label = "Mean")
    plot!(title = "Non-linear signal", xlab = "t [sec]", ylab = "η [m]")

    pltb = plot(Ltₑᵥ, Levents, legend = false, grid=true)
    plot!(tmeanLNR,LNR_mean,line = :dash, color = :red, linewidth = 2, label = "Mean")
    plot!(title = "Linear signal", xlab = "t [sec]", ylab = "η [m]")

    plt_ab = plot(plta, pltb, layout = @layout [a ; b ])
end

end
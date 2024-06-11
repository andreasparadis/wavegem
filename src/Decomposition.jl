module Decomposition
using Plots, LaTeXStrings, DelimitedFiles, Dates, LinearAlgebra
using FFTW, DSP, Statistics, BSplineKit

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
#############################################################################################
# Make directories
make_dirs(2, Decpath, DecFigs)
make_dirs(2, OFASTpath*"/ExtElev/$(case_id)", OFASTpath*"/ExtElev/$(case_id)/$(run_id)")

#############################################################################################
# Import time-domain data for each pi/2 phase shifted realization
## Reference signal
fid = Decpath*sigf[1] # Path to file
open(fid, "r")
A00 = readdlm(fid, '\t', Float64, '\n')

## Phase shifted signals
∅ = zeros(Float64, length(A00[:,1]),3);     A = ∅[:,:]
for i ∈ 1:3
    fid = Decpath*sigf[i+1] # Path to file
    open(fid, "r")
    cont = readdlm(fid, '\t', Float64, '\n')
    A[:,i] = cont[:,2] 
end
A05 = A[:,1];   A10 = A[:,2];   A15 = A[:,3] 

#############################################################################################
# Signal truncation
## Truncate time vector
dt = A00[2,1]-A00[1,1]                  # Time step
ibeg = Int64(round(tstart/dt))+1        # lower limit
iend = Int64(round(tend/dt))            # upper limit
t00 = A00[ibeg:iend, 1] .- A00[ibeg, 1] # Set first value as time zero
tOG = round.(t00*100)/100               # Global time vector for all components
nₜ = length(tOG)

## Truncate vectors of values
A̅00 = A00[ibeg:iend, 2]
A̅05 = A05[ibeg:iend]
A̅10 = A10[ibeg:iend]
A̅15 = A15[ibeg:iend]

#############################################################################################
# 3rd order decomposition of original signal (signal 0)
SPEC0, AMPS, ARGS, sig_comps, sig_recon, L, fₚ = decomp_3rd(tOG, A̅00, A̅05, A̅10, A̅15)
FR1, Â₀₀, ϕ₀₀ = SPEC0[:,1], SPEC0[:,2], SPEC0[:,3] # Spectral info of signal 0

#############################################################################################
# Isolate, shift around 0 and plot 'pulses' based on peaks
## Find peaks and sort them in descending order
## Isolate pulses given a specified time range (increase the range to consider more pusles)
t_range = 50  # [s]
NP = 8
shift = 2
res_f = 100

MinPeakVal = 2*std(A̅00)
MinPeakDist = WAVEGEM.Tᵢ
tₑᵥⁱ, ev_int, tₑᵥ, events, rng, Tₛᵉᵛ, PeakId, PeakPos, PeakVal = 
    ev_identify(A̅00, tOG, MinPeakVal, MinPeakDist, t_range, NP, shift, res_f)

Val_LNR = sig_comps[:,1]
MinPeakVal = 2*std(Val_LNR)
Ltₑᵥⁱ, Lev_int, Ltₑᵥ, Levents, rng, LTₛᵉᵛ, L_PeakId, L_PeakPos, L_PeakVal = 
    ev_identify(Val_LNR, tOG, MinPeakVal, MinPeakDist, t_range, NP, shift, res_f)

## Calculate mean pulse
timemean = tₑᵥ[:,NP]
NLN_mean = mean!(ones(2*rng+1,1), events)

tmeanLNR = Ltₑᵥ[:,NP]
LNR_mean = mean!(ones(2*rng+1,1), Levents)

#############################################################################################
# Compare pulses and calculate statistical measures
t_rng = 5  # [s]
id_rng = round(Int, (t_rng/2) / dt + 1)  # +- in terms of time intervals
ID0 = div(size(Levents, 1), 2)  # Index of t=0 (middle of vector)

## Pulses corresponding to the shortened time range
shrt = Levents[ID0-id_rng:ID0+id_rng, :]

# Find the pulse with the highest peak
M_pls_i = findmax(shrt)[2]
fmf, inx = M_pls_i[1], M_pls_i[2]

## 'inx' is the identifier of the pulse containing the overall max
A = shrt[:, inx]

## Statistical comparison of other pulses against the identified pulse
## Exclude pulse A from the set of pulses
B = hcat(shrt[:, 1:inx-1], shrt[:, inx+1:end])

cvar = cov(A, B)
sigA = std(A)
sigB = std(B, dims=1)

## Pearson correlation coefficient
rho = cvar ./ (sigA * sigB)

R = A .- B
Ravg = mean(R, dims=1)
sigR = std(R, dims=1)

#############################################################################################
# Examine most extreme event
lb = PeakId[1]-500; ub = lb+1000 # For event plots

Aₑᵥ = ev_int[:,1]
tev = tₑᵥⁱ[:,1]

freq, mag, _,_,_,Nfft, H = one_side_asp(Aₑᵥ, tev)
itp_sp = interpolate(freq, mag, BSplineOrder(4))
FR = range(0,fcut, Nfft)
MAG = itp_sp.(FR)
EV = ifft(H)

plt = plot(freq, mag)
plot!(FR, MAG)
plot!(xlim=(0,fcut))
display(plt)
#############################################################################################
if frec
    # Surface elevation (temporal)
    fid_elev::String = "eta_lin" # File name
    open(Decpath*fid_elev, "w")
    writedlm(Decpath*fid_elev, [tOG Val_LNR], '\t')

    # Surface elevation spectrum
    fid_elev::String = "eta_lin_spec" # File name
    open(Decpath*fid_elev, "w")
    writedlm(Decpath*fid_elev, [FR1 AMPS[:,1]], '\t')

    # # Event
    # evdir::String = "events/"
    # # if isdir(evdir)
    # #     error("Andreas: This folder already exists, since this event has already been recorded!")
    # # else
    # #     mkdir("library/"*Decpath*evdir)
    # # end

    # fid_ev::String = "event_1" # File name
    # open(Decpath*evdir*fid_ev, "w")
    # writedlm(Decpath*evdir*fid_ev, [tOG[lb:ub] A̅00[lb:ub]], '\t')

    # # lbe = Int(round(360/dt)); ube = Int(round(415/dt))
    # # fid_ev::String = "event_2" # File name
    # # open("library/"*Decpath*evdir*fid_ev, "w")
    # # writedlm("library/"*Decpath*evdir*fid_ev, [tOG[lbe:ube] A̅00[lbe:ube]], '\t')

    # fid_ev::String = "event_1_int" # File name
    # open(Decpath*evdir*fid_ev, "w")
    # writedlm(Decpath*evdir*fid_ev, [tOG Aₑᵥ], '\t')
end
#############################################################################################
## PLOTS 
if fplot
    # For JONSWAP spectrum comparisons
    Tₑ₀ = sqrt(2π/g * Ldom/2) # Initial guess for maximum period
    long_wave = wave_qnts(Tₑ₀,d)
    ω⁻ = long_wave.ω
    Tₑ = round(2π/ω⁻) # Maximum period

    Hₛ = 4*std(A̅00)
    ω⁺ = 2π/WAVEGEM.Tᵢ
    dω = 2π/((nₜ-1)*dt)
    Ncut = Int64(round((ω⁺-ω⁻)/dω))
    Tₚ = 1/fₚ

    fⱼ,Sⱼ = spectrum(Hₛ,Tₚ,γ,WAVEGEM.Tᵢ,Tₑ,Ncut) # Generate JONSWAP spectrum

    # Plot decomposed signal and its components
    plt1 = plot(tOG[lb:ub], Val_LNR[lb:ub], xlab = "Time [sec]", ylab = "Elevation [m]", line = :solid, color = :blue, label = "Linear")
    plot!(tOG[lb:ub], sig_comps[lb:ub,2], line = :solid, color = :green, label = "2nd-")
    plot!(tOG[lb:ub], sig_comps[lb:ub,3], line = :solid, color = :magenta, label = "2nd+")
    plot!(tOG[lb:ub], sig_comps[lb:ub,4], line = :solid, color = :black, label = "3rd")
    plot!(tOG[lb:ub], A̅00[lb:ub], line = :dash, color = :red, linewidth = 1, label = "Original")
    display(plt1)
    savefig(DecFigs*"decomp_elev.svg")
    savefig(DecFigs*"decomp_elev.png")

    plt2 = plot(tOG[lb:ub], [sig_recon[lb:ub] A̅00[lb:ub]], xlab = "Time [sec]", ylab = "Elevation [m]", line = [:solid :dash], color = [:blue :red], label = ["Sum" "Original"])
    savefig(DecFigs*"elev_compare.svg")
    savefig(DecFigs*"elev_compare.png")

    # plt3 = plot(tₑᵥ, Aₑᵥ, line = :solid, color = :blue, linewidth = 1, label = "Interpolation")
    # plot!(tOG[lb:ub], A̅00[lb:ub], line = :dash, color = :red, linewidth = 1, label = "Original")
    # display(plt3)

    # Plot identified peaks on top original and linear signals
    plot(tOG, A̅00, xlab = "Time [sec]", ylab = "Elevation [m]",  label = "Original signal")
    display(plot!(PeakPos, PeakVal, seriestype=:scatter, ms=2, mc=:red, lab = "Peaks",legend=:outerbottom, legendcolumns=2))
    savefig(DecFigs*"elev_peaks_nln.svg")
    savefig(DecFigs*"elev_peaks_nln.png")

    plot(tOG, Val_LNR, xlab = "Time [sec]", ylab = "Elevation [m]",  label = "Linear signal")
    display(plot!(L_PeakPos, L_PeakVal, seriestype=:scatter, ms=2, mc=:red, lab = "Peaks",legend=:outerbottom, legendcolumns=2))
    savefig(DecFigs*"elev_peaks_lin.svg")
    savefig(DecFigs*"elev_peaks_lin.png")

    # Plot 'pulses' and mean
    plot(tₑᵥ, events, title="Nonlinear Signal", xlabel="Time [sec]", ylabel="Elevation [m]", grid=true)
    plta = plot!(timemean,NLN_mean,line = :dash, color = :red, linewidth = 2, label = "Mean")

    plot(Ltₑᵥ, Levents, title="Linearized Signal", xlabel="Time [sec]", ylabel="Elevation [m]", grid=true)
    pltb = plot!(tmeanLNR,LNR_mean,line = :dash, color = :red, linewidth = 2, label = "Mean")

    display(plot(plta, pltb, layout = @layout [a ; b ]))
    savefig(DecFigs*"elev_pulses_compare.svg")
    savefig(DecFigs*"elev_pulses_compare.png")

    # Plot spectral decomposition
    df = (FR1[end]-FR1[1])/(L/2);
    dfⱼ = abs(fⱼ[end]-fⱼ[1])/(length(fⱼ)-1)
    η̂ = sqrt.(2*Sⱼ*dfⱼ)

    flb = 1; fub = Int(round(fcut/df))
    plt_sp = plot(FR1[flb:fub], Â₀₀[flb:fub], xlab = "f [Hz]", ylab = L"S(f)~[m^2 s]", line =:dashdot, lab="Original", yscale=:log10)
    plot!(FR1[flb:fub], AMPS[flb:fub,1], line =:solid, lw=0.5, label = "Linear", yscale=:log10)
    plot!(FR1[flb:fub], AMPS[flb:fub,2], line =:solid, lw=0.5, label = "2nd-", yscale=:log10)
    plot!(FR1[flb:fub], AMPS[flb:fub,3], line =:dash, lw=0.5, label = "2nd+", yscale=:log10)
    plot!(FR1[flb:fub], AMPS[flb:fub,4], line =:dot, lw=0.5, label = "3rd", yscale=:log10)
    plot!(fⱼ, η̂, lab="JONSWAP", color=:red, lw=2, yscale=:log10)
    plot!(xlim=(0,fcut), ylim=(1e-6,1), minorgrid=true)
    display(plt_sp)
    savefig(DecFigs*"decomp_spectra.svg")
    savefig(DecFigs*"decomp_spectra.png")

    freq, mag, _ = one_side_asp(Aₑᵥ,tev)
    plt_ev = plot(FR1[flb:fub], Â₀₀[flb:fub], xlab = "f [Hz]", ylab = L"S(f)~[m^2 s]", line =:dashdot, lab="Original", yscale=:log10)
    plot!(freq, mag*res_f, line =:solid, lab="Event", yscale=:log10)
    plot!(fⱼ, η̂, lab="JONSWAP", color=:red, lw=2, yscale=:log10)
    plot!(xlim=(0,fcut), ylim=(1e-6,1), minorgrid=true)
    display(plt_ev)
    savefig("library/"*Decpath*evdir*"event_spectra.svg")
    savefig("library/"*Decpath*evdir*"event_spectra.png")

    # Compare spectra of input signal and OW3D simulation
    elev_fid::String = "library/"*rundir*"0/eta_t"
    open(elev_fid, "r")   
    elev = readdlm(elev_fid, '\t', Float64, skipstart=1, '\n')
    tinp, ηinp = elev[:,1], elev[:,2]
    freq, mag, _ = one_side_asp(ηinp,tinp)

    plt_comp = plot(freq, mag, xlab = "f [Hz]", ylab = L"S(f)~[m^2 s]", lab="Input", yscale=:log10)
    plot!(FR1[flb:fub], Â₀₀[flb:fub], line =:dot, lab="Original", yscale=:log10)
    plot!(FR1[flb:fub], AMPS[flb:fub,1]/L, line =:dashdot, lw=0.2, label = "Linear", yscale=:log10)
    plot!(fⱼ, η̂, lab="JONSWAP", color=:red, lw=2, yscale=:log10)
    plot!(xlim=(0,fcut), ylim=(1e-6,1), minorgrid=true)
    display(plt_comp)
    savefig(DecFigs*"spect_comp.svg")
    savefig(DecFigs*"spect_comp.png")

    h = PeakVal[1]
    plt_groups = plot(tOG, A̅00, xlab = "Time [sec]", ylab = "Elevation [m]",  label = "Original signal")
    for i ∈ 1:max(20, length(PeakVal))
        c = PeakPos[i]
        w = 25
        plot!(Shape([c-w,c-w,c+w,c+w],[-h,h,h,-h]), color=:red, opacity=0.2, label = "WG$i")
    end
    plot!(legend = false)
    display(plt_groups)
    savefig(DecFigs*"WGs.svg")

    plt_groups = plot(tOG, A̅00, xlab = "Time [sec]", ylab = "Elevation [m]",  label = "Original signal")

    c = PeakPos[1]
    w = 50
    plot!(Shape([c-w,c-w,c+w,c+w],[-h,h,h,-h]), color=:red, opacity=0.3, label = "Event")
    plot!(xlim=(6000,7000))
    display(plt_groups)
    savefig(DecFigs*"Event1.svg")
end
# # Find zero Up or Down Crossings
# UCiD = findall(diff(sign.(A̅00)) .== 2)  # Up-crossing (-2 for Down-crossing)
# UC_Val = A̅00[UCiD]
# UC_Pos = tOG[UCiD]

# UCiD = findall(diff(sign.(Val_LNR)) .== 2)  # Up-crossing (-2 for Down-crossing)
# L_UC_Val = Val_LNR[UCiD]
# L_UC_Pos = tOG[UCiD]

# # Plot original signal along with up crossing points
# plt = plot(tOG, A̅00, line=:solid, color=:blue, lab="Elevation")
# plot!(UC_Pos, UC_Val, seriestype=:scatter, mc=:red, ms=2, lab = "Up-Crossing points")
# plot!(
#     xlim=(tstart, tstart+200),
#     xlabel="Time [sec]", ylabel="Elevation [m]",
#     legend=:topright
# )
# display(plt)
end
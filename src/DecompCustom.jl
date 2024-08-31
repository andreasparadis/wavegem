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
## Significant wave height [m], Peak period [s], peakedness [-], Cut-off frequency [Hz]
Hₛ, Tₚ, γ, fcut::Float64 = 0.069, 1, 3.3, 4  # JONSWAP parameters 
const ρ, g::Float64 = 1025.0, 9.81  # Sea water density [kg/m³], Accel. of gravity [m/s²]
Ldom, d::Float64 = 10.0, 0.4    # [m] Domain length, Water depth

prb_id::Int8 = 4    # Column No corresponding to probe at phase focusing location


# Decpath =  "/home/andreasp/WAVEGEM/library/JFM/EV_2/G0/CB_DTA/"
Decpath =  "/home/andreasp/WAVEGEM/library/JFM/FULL/"
DecFigs = Decpath

tstart = 0; tend = 512
frec, fplot = Bool(1), Bool(1)
# sigf = ("A00.txt", "A05.txt", "A10.txt", "A15.txt") # Signal files for decomposition
sigf = ("FR1_1_0.txt", "FR1_1_piover2.txt", "FR1_1_pi.txt", "FR1_1_3piover2.txt") # Signal files for decomposition

#############################################################################################
# Import t-domain data for each pi/2 phase shifted realization
## Reference signal
fid = Decpath*sigf[1] # Path to file
fcont = parse_fxw(fid, 0)
A00 = [fcont[:,1] fcont[:,prb_id]]

## Phase shifted signals
∅ = zeros(Float64, length(A00[:,1]),3);     A = ∅[:,:]
for i ∈ 1:3
    fid = Decpath*sigf[i+1] # Path to file
    cont = parse_fxw(fid, 0)
    A[:,i] = cont[:,prb_id] 
end
A05 = A[:,1];   A10 = A[:,2];   A15 = A[:,3] 

#############################################################################################
# Signal truncation
## Truncate t vector
dt = A00[2,1]-A00[1,1]                  # t step
ibeg = Int64(round(tstart/dt))+1        # lower limit
iend = Int64(round(tend/dt))            # upper limit
t00 = A00[ibeg:iend, 1] .- A00[ibeg, 1] # Set first value as t zero
tOG = round.(t00*100)/100               # Global t vector for all components
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
df = (FR1[end]-FR1[1])/(L/2)

#############################################################################################
# Isolate, shift around 0 and plot 'pulses' based on peaks
## Find peaks and sort them in descending order
## Isolate peak sections given a specified t range
t_range = 20  # [s]
iNP = 8
shift = 2
res_f = 100

MinPeakVal = 2*std(A̅00)
MinPeakDist = 1/fcut

tₑᵥⁱ, ev_int, tₑᵥ, events, rng, Tₛᵉᵛ, PeakId, PeakPos, PeakVal, NP = 
    ev_identify(A̅00, tOG, MinPeakVal, MinPeakDist, t_range, iNP, shift, res_f)

Val_LNR = sig_comps[:,1]
MinPeakVal = 2*std(Val_LNR)

Ltₑᵥⁱ, Lev_int, Ltₑᵥ, Levents, rng, LTₛᵉᵛ, L_PeakId, L_PeakPos, L_PeakVal, LNP = 
    ev_identify(Val_LNR, tOG, MinPeakVal, MinPeakDist, t_range, iNP, shift, res_f)

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
# Examine most extreme event
Aₑᵥ = ev_int[:,1]
tev = tₑᵥⁱ[:,1]

freq, mag, _,_,_,Nfft, H = one_side_asp(Aₑᵥ, tev)
itp_sp = interpolate(freq, mag, BSplineOrder(4))
FR = range(0,fcut, Nfft)
MAG = itp_sp.(FR)
EV = ifft(H)

plt_spitp = plot(freq, mag, lw=2, xlab = "f [Hz]", ylab = L"Amplitude", lab="Max event")
plot!(FR, MAG, line =:dashdot, lw=2, lab="Interpolation")
plot!(xlim=(0,fcut))

#############################################################################################
if frec
    # Surface η (temporal)
    fid_elev::String = "eta_lin" # File name
    open(Decpath*fid_elev, "w")
    writedlm(Decpath*fid_elev, [tOG Val_LNR], '\t')

    # Surface η spectrum
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
    # writedlm(Decpath*evdir*fid_ev, [tpart A̅00[lb:ub]], '\t')

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
    ω⁺ = 2π*fcut
    dω = 2π/((nₜ-1)*dt)
    Ncut = Int64(round((ω⁺-ω⁻)/dω))
    Tₚ = 1/fₚ

    fⱼ,Sⱼ = spectrum(Hₛ,Tₚ,γ,1/fcut,Tₑ,Ncut) # Generate JONSWAP spectrum
    dfⱼ = abs(fⱼ[end]-fⱼ[1])/(length(fⱼ)-1)
    η̂ = sqrt.(2*Sⱼ*dfⱼ)

    # For partial plots
    lb = PeakId[1]-rng; ub = lb+2*rng 
    tpart = tOG[lb:ub]

    # 1: Decomposed signal and its components
    plt1 = plot(tpart, A̅00[lb:ub], lw = 2, lab = "Non-linear")
    plot!(tpart, Val_LNR[lb:ub], lw = 2, label = L"\mathbf{1_{st}}")
    plot!(tpart, sig_comps[lb:ub,2], line =:dot, lw = 2, lab = L"\mathbf{2_{nd}^{-}}")
    plot!(tpart, sig_comps[lb:ub,3], line =:dot, lw = 2, lab = L"\mathbf{2_{nd}^{+}}")
    plot!(tpart, sig_comps[lb:ub,4], line =:dot, lw = 2, lab = L"\mathbf{3_{rd}}")
    plot!(title = "Decomposed signal components", xlab = "t [sec]", ylab = "η [m]")
    # plot!(legendcolumns=3)

    # 2: Reconstructed signal against the original
    plt_recon = plot(tpart, A̅00[lb:ub], line =:solid, label = "Non-linear")
    plot!(tpart, sig_recon[lb:ub], line =:dash, color =:red, label = L"\sum components")
    plot!(title = "Comparison: Original vs Reconstructed Signal", xlab = "t [sec]", ylab = "η [m]")

    # 3: Interpolated event against original event
    plt_itp = plot(tev, Aₑᵥ, line = :solid, lw = 1, lab = "Interpolation")
    plot!(twiny(), tpart, A̅00[lb:ub], line =:dash, color =:red, lw = 1, ylab = "η [m]", lab = "Non-linear")
    plot!(title = "Event at Highest Peak: Original vs Interpolated", xlab = "t [sec]")

    # 4: Peaks on top of original signal
    plt_peaks = plot(tOG, A̅00, label = "Non-linear")
    plot!(PeakPos, PeakVal, seriestype=:scatter, ms=2, mc=:red, lab = "Peaks",legend=:bottomleft, legendcolumns=2)
    plot!(title = "Peaks of Non-linear signal", xlab = "t [sec]", ylab = "η [m]")

    # 5: Peaks on top of linear signal
    plt_Lpeaks = plot(tOG, Val_LNR, xlab = "t [sec]", ylab = "η [m]",  label = "Linear")
    plot!(L_PeakPos, L_PeakVal, seriestype=:scatter, ms=2, mc=:red, lab = "Peaks",legend=:bottomleft, legendcolumns=2)
    plot!(title = "Peaks of Linear signal", xlab = "t [sec]", ylab = "η [m]")

    # 6: Peak sections and mean of sections
    plta =  plot(tₑᵥ, events, legend = false, grid=true)
    plot!(tmean,NLN_mean,line = :dash, color = :red, linewidth = 2, label = "Mean")
    plot!(title = "Non-linear signal", xlab = "t [sec]", ylab = "η [m]")

    pltb = plot(Ltₑᵥ, Levents, legend = false, grid=true)
    plot!(tmeanLNR,LNR_mean,line = :dash, color = :red, linewidth = 2, label = "Mean")
    plot!(title = "Linear signal", xlab = "t [sec]", ylab = "η [m]")

    plt_ab = plot(plta, pltb, layout = @layout [a ; b ])

    # 7: Spectral decomposition
    flb = 1; fub = Int(round(fcut/df))

    plt_sp = plot(FR1[flb:fub], Â₀₀[flb:fub], line =:dashdot, lab="Non-linear")
    plot!(FR1[flb:fub], AMPS[flb:fub,1], line =:solid, lw=0.5, label = L"\mathbf{1_{st}}")
    plot!(FR1[flb:fub], AMPS[flb:fub,2], line =:solid, lw=0.5, label = L"\mathbf{2_{nd}^{-}}")
    plot!(FR1[flb:fub], AMPS[flb:fub,3], line =:dash, lw=0.5, label = L"\mathbf{2_{nd}^{+}}")
    plot!(FR1[flb:fub], AMPS[flb:fub,4], line =:dot, lw=0.5, label = L"\mathbf{3_{rd}}")
    plot!(fⱼ, η̂, lab="JONSWAP", color=:red, lw=2)
    plot!(title = "Spectral Decomposition", xlab = L"f~[Hz]", ylab = L"S(f)~[m^2 s]")
    plot!(yscale=:log10, legendcolumns=1)
    plot!(xlim=(0,fcut), ylim=(1e-6,1), minorgrid=true)

    # 8:
    freq, mag, _ = one_side_asp(Aₑᵥ,tev)
    plt_evsp = plot(FR1[flb:fub], Â₀₀[flb:fub], line =:dashdot, lab="Non-linear")
    plot!(freq, mag, line =:solid, lab="Event")
    plot!(fⱼ, η̂, lab="JONSWAP", color=:red, lw=2)
    plot!(yscale=:log10, legendcolumns=1)
    plot!(xlim=(0,fcut), ylim=(1e-6,1), minorgrid=true)
    plot!(title = "Comparison of Spectra", xlab = L"f~[Hz]", ylab = L"S(f)~[m^2 s]")
   
    # 9: 
    h = PeakVal[1]
    plt_groups = plot(tOG, A̅00, xlab = "t [sec]", ylab = "η [m]",  label = "Non-linear")
    for i ∈ 1:min(20, length(PeakVal))
        c = PeakPos[i]
        w = 25
        plot!(Shape([c-w,c-w,c+w,c+w],[-h,h,h,-h]), color=:red, opacity=0.2, label = "WG$i")
    end
    plot!(legend = false)
    
    # 11:
    c = PeakPos[1]; w = 50
    plt_maxev = plot(tOG, A̅00, xlab = "t [sec]", ylab = "η [m]",  label = "Non-linear")
    plot!(Shape([c-w,c-w,c+w,c+w],[-h,h,h,-h]), color=:red, opacity=0.3, label = "Event")
    plot!(xlim=(c-10*w, c+10*w))
    
    # Plots to display
    display(plt1)
    display(plt_recon)
    display(plt_itp)
    display(plt_peaks)
    display(plt_Lpeaks)
    display(plt_ab)
    display(plt_sp)
    display(plt_evsp)
    display(plt_groups)
    display(plt_maxev)
    display(plt_spitp)

    # Figures to save
    savefig(plt1, DecFigs*"decomp_elev.png");   savefig(plt1, DecFigs*"decomp_elev.svg")
    savefig(plt_recon, DecFigs*"elev_compare.png"); savefig(plt_recon, DecFigs*"elev_compare.svg")
    savefig(plt_peaks, DecFigs*"elev_peaks_nln.png"); savefig(plt_peaks, DecFigs*"elev_peaks_nln.svg")
    savefig(plt_Lpeaks, DecFigs*"elev_peaks_lin.png"); savefig(plt_Lpeaks, DecFigs*"elev_peaks_lin.svg")
    savefig(plt_ab, DecFigs*"elev_pulses_compare.png"); savefig(plt_ab, DecFigs*"elev_pulses_compare.svg")
    savefig(plt_sp, DecFigs*"decomp_spectra.png"); savefig(plt_sp, DecFigs*"decomp_spectra.svg")
    # savefig(plt_evsp, Decpath*evdir*"event_spectra.png"); savefig(plt_evsp, Decpath*evdir*"event_spectra.svg")
    savefig(plt_groups, DecFigs*"WGs.svg")
    savefig(plt_maxev, DecFigs*"Event1.svg")
end
# # Find zero Up or Down Crossings
# UCiD = findall(diff(sign.(A̅00)) .== 2)  # Up-crossing (-2 for Down-crossing)
# UC_Val = A̅00[UCiD]
# UC_Pos = tOG[UCiD]

# UCiD = findall(diff(sign.(Val_LNR)) .== 2)  # Up-crossing (-2 for Down-crossing)
# L_UC_Val = Val_LNR[UCiD]
# L_UC_Pos = tOG[UCiD]

# # Plot original signal along with up crossing points
# plt = plot(tOG, A̅00, line=:solid, color=:blue, lab="η")
# plot!(UC_Pos, UC_Val, seriestype=:scatter, mc=:red, ms=2, lab = "Up-Crossing points")
# plot!(
#     xlim=(tstart, tstart+200),
#     xlabel="t [sec]", ylabel="η [m]",
#     legend=:topright
# )
# display(plt)
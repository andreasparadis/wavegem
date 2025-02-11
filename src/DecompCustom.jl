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
prb_id::Int8 = 4        # Probe id at phase focusing location (Channel No)
prb_pos::Float64 = 8.0  # Probe position [m]
tdur::Float64 = 1200.0  # Signal time range to consider [s]
evID::Int8 = 5          # id of event to save (desc. order from heighest peak)

# Data locations - Edit accordingly
srcfid = joinpath("HS006","ETA_HS006TP1TR128.txt")  # Theoretical time series
cfld::String = joinpath("HS006TP1TR128","Random")   # Case folder (Experimental)
# cfld = joinpath("HS006TP1TR128","Random","CORR")  # Case folder (Corrected Exp.)

sigf = ("A00.txt", "A05.txt", "A10.txt", "A15.txt") # Signals to decompose

# Paths (no need to edit)
libpath = joinpath(pwd(),"library","UCL","Waves") # Library folder
srcpath = joinpath(libpath,"AP_FR",srcfid)
Parent = joinpath(libpath,cfld)
Decpath = joinpath(Parent,"Decomposition")
CBpath = joinpath(Parent,"CLB")
DecFigs = joinpath(Decpath,"Figures")
# evdir = joinpath(Decpath,"EV$evID")
evdir = joinpath(Decpath,"Surge_EV$evID")

# Flags
fevid = true  # Event identification (true only for the complete run)
fplot = true   # Plot results
frec = true    # Write output files & figures
tshift = true  # True for experimental data

#############################################################################################
# Calibrate data the first time a new set of runs is examined
if !isdir(CBpath)
    for i ∈ 1:4
        cb_exp(Parent,i) # (Folder with exp. data, phase shift id)
    end
end
# Make other necessary directories
if !isdir(Decpath)
    mkdir(Decpath)
end
if !isdir(DecFigs)
    mkdir(DecFigs)
end
if !isdir(evdir)
    mkdir(evdir)
end
#############################################################################################
if !@isdefined fcont # To not re-read the txt inputs if already loaded
    # Import t-domain data for each pi/2 phase shifted realization
    ## Reference signal
    fid = joinpath(CBpath,sigf[1]) # Path to file
    fcont = parse_fxw(fid, 1)
    A00 = [fcont[:,1] fcont[:,prb_id+1]]

    ## Phase shifted signals
    ∅ = zeros(Float64, length(A00[:,1]),3);     A = ∅[:,:]
    for i ∈ 1:3
        fid2 = joinpath(CBpath,sigf[i+1])
        cont = parse_fxw(fid2, 1)
        A[:,i] = cont[:,prb_id+1] 
    end
    A05 = A[:,1];   A10 = A[:,2];   A15 = A[:,3] 

    ## Theoretical signal used in the wavemaker
    srccont = parse_fxw(srcpath, 3)
    tˢ = srccont[:,1]
    ηˢ = srccont[:,2]

    FRˢ, MAGˢ, ϕˢ,_ = one_side_asp(ηˢ, tˢ)
end

##########################################################################################
# Signal truncation
## Truncate t vector
Tᵢ = 1/fcut             # Cut-off period of spectrum [s]
λ⁻ = g/(2π) * Tᵢ^2     # Shortest wavelength (Deep water) [m]
υ⁻ = λ⁻/Tᵢ              # Slowest wave [m/s]
tstart = prb_pos/υ⁻; tend = tstart+tdur

dt = A00[2,1]-A00[1,1]                  # t step
ibeg = Int64(round(tstart/dt))+1        # lower limit
iend = Int64(round(tend/dt))            # upper limit
t00 = A00[ibeg:iend, 1] .- A00[ibeg, 1] # Set first value as t zero
tOG = round.(t00*100)/100               # Global t vector for all components
nₜ = length(tOG)

## Truncate vectors of values
A̅00 = A00[ibeg:iend, 2] #.- mean(A00[:,2])
A̅05 = A05[ibeg:iend] #.- mean(A05)
A̅10 = A10[ibeg:iend] #.- mean(A10)
A̅15 = A15[ibeg:iend] #.- mean(A15)

#############################################################################################
# 3rd order decomposition of original signal (signal 0)
SPEC0, AMPS, ARG, sig_comps, sig_recon, L, fₚ, pltϕ = decomp_3rd(tOG, A̅00, A̅05, A̅10, A̅15)
# SPEC0, AMPS, ARG, sig_comps, sig_recon, L, fₚ = decomp_3rd(tOG, A̅00, A̅05, A̅10, A̅15)

FR1, Â₀₀, ϕ₀₀ = SPEC0[:,1], SPEC0[:,2], SPEC0[:,3] # Spectral info of signal 0
df = (FR1[end]-FR1[1])/(L/2)
Val_LNR = sig_comps[:,1]    # 1st order component

# Instantaneous frequency
𝓗 = hilbert(Val_LNR)     
u = abs.(𝓗)
θ = angle.(𝓗);  θ = unwrap(θ)

omH = zeros(Float64, nₜ)
tom = zeros(Float64, nₜ)
for i ∈ 2:nₜ-1
    omH[i] = (θ[i+1] - θ[i-1]) / (tOG[i+1]-tOG[i-1])
    if omH[i] > 0 
        tom[i] = 2π
    else
        tom[i] = 0
    end
end
#############################################################################################
if fevid
    # Event Identification
    ## Find peaks and sort them in descending order
    ## Isolate peak sections given a specified t range
    t_range = 10  # [s]
    iNP = 8
    shift = 2
    res_f = 10

    cstd = 3
    MinPeakDist = 1/fcut
    MinPeakVal = cstd*std(A̅00)

    tₑᵥⁱ, ev_int, tₑᵥ, events, rng, Tₛᵉᵛ, PeakId, PeakPos, PeakVal, NP = 
        ev_identify(A̅00, tOG, MinPeakVal, MinPeakDist, t_range, iNP, shift, res_f)

    LMinPeakVal = cstd*std(Val_LNR)

    Ltₑᵥⁱ, Lev_int, Ltₑᵥ, Levents, rng, LTₛᵉᵛ, L_PeakId, L_PeakPos, L_PeakVal, LNP = 
        ev_identify(Val_LNR, tOG, LMinPeakVal, MinPeakDist, t_range, iNP, shift, res_f)

    ###################################################################################
    # Select event duration using instantaneous frequency
    ## Event limits
    lb = PeakId[evID]
    while sign(omH[lb]) == sign(omH[lb-1])
        global lb = lb - 1
    end

    while sign(Val_LNR[lb-1])-sign(Val_LNR[lb]) ≠ 2
        global lb = lb - 1
    end

    ub = PeakId[evID]
    while sign(omH[ub]) == sign(omH[ub+1])
        global ub = ub + 1
    end

    while sign(Val_LNR[ub+1])-sign(Val_LNR[ub]) ≠ 2
        global ub = ub + 1
    end

    # Examine extreme event
    # lb = PeakId[evID]-rng; ub = lb+2*rng 
    tst = 1055;	ten = 1075
    lb = Int64(round(tst/dt))
    ub = Int64(round(ten/dt))

    tpart = tOG[lb:ub]; Lₚ = length(tpart)

    Aₑₑ = A̅00[lb:ub]
    Aₑₑ¹ = Val_LNR[lb:ub]
    tₑₑ = tOG[lb:ub] .- tOG[lb]

    # Aₑₑ = ev_int[:,evID]
    # tₑₑ = tₑᵥⁱ[:,evID]

    fr_ev, mag_ev, phi,_,_,Nfft, H = one_side_asp(Aₑₑ¹, tₑₑ)
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
    # 'Pulses' (shorter section around event peak)
    ## Calculate mean pulse
    tmean = tₑᵥ[:,NP]
    NLN_mean = mean!(ones(2*rng+1,1), events)

    tmeanLNR = Ltₑᵥ[:,LNP]
    LNR_mean = mean!(ones(2*rng+1,1), Levents)

    # Compare signal section and calculate statistical measures
    sslen = 5  # Signal section length[s]
    ρᴾ, R̅, σᴿ = sect_corr(sslen, dt, Levents)
else
    lb = 1
    ub = length(Val_LNR)
    tpart = tOG[lb:ub]; Lₚ = length(tpart)
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

    fⱼ,Sⱼ = spectrum(Hₛ,Tₚ,γ,Tᵢ,Tₑ,Ncut) # Generate JONSWAP spectrum
    dfⱼ = abs(fⱼ[end]-fⱼ[1])/(length(fⱼ)-1)
    η̂ = sqrt.(2*Sⱼ*dfⱼ)

    # 1: Spectral decomposition
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

    # 2: Decomposed signal and its components
    plt1 = plot(tpart, A̅00[lb:ub], lw = 2, lab = "Non-linear")
    plot!(tpart, Val_LNR[lb:ub], lw = 2, label = L"\mathbf{1_{st}}")
    plot!(tpart, sig_comps[lb:ub,2], line =:dot, lw = 2, lab = L"\mathbf{2_{nd}^{-}}")
    plot!(tpart, sig_comps[lb:ub,3], line =:dot, lw = 2, lab = L"\mathbf{2_{nd}^{+}}")
    plot!(tpart, sig_comps[lb:ub,4], line =:dot, lw = 2, lab = L"\mathbf{3_{rd}}")
    plot!(title = "Decomposed signal components", xlab = "t [sec]", ylab = "η [m]")
    # plot!(legendcolumns=3)

    # 3: Reconstructed signal against the original
    plt_recon = plot(tpart, A̅00[lb:ub], line =:solid, label = "Non-linear")
    plot!(tpart, sig_recon[lb:ub], line =:dash, color =:red, label = L"\sum components")
    plot!(title = "Comparison: Original vs Reconstructed Signal", xlab = "t [sec]", ylab = "η [m]")

    # 4: Compare theoretical and resulting experimental time series
    plt_tcomp = plot(xlab = "t [sec]", ylab = "η [m]")
    plot!(tˢ, ηˢ,  lab = "Theoretical elevation")
    plot!(tOG, A̅00, lab = "Experimental (probe $prb_id)")

    # 5: Compare theoretical and resulting experimental spectra
    plt_spcomp = plot(xlab = "f [Hz]", ylab = "S [m]", xlim=(0,fcut))
    plot!(FRˢ, MAGˢ,  lab = "Theoretical amp. spectrum")
    plot!(FR1, Â₀₀, lab = "Experimental (probe $prb_id)",line=:dot,opacity=0.7)
    plot!(FR1, AMPS[:,1], lab="1st order",line=:dot,opacity=0.5)

    if fevid
        # 6: Peaks on top of original signal
        plt_peaks = plot(tOG, A̅00, label = "Non-linear")
        plot!(PeakPos[1:NP], PeakVal[1:NP], seriestype=:scatter, ms=2, mc=:red, lab = "Peaks",legend=:bottomleft, legendcolumns=3)
        plot!([tOG[1];tOG[end]],[MinPeakVal;MinPeakVal], lab="$cstd η", mc=:red, ls=:dash)
        plot!(title = "Peaks of Non-linear signal", xlab = "t [sec]", ylab = "η [m]")

        # 7: Peaks on top of linear signal
        plt_Lpeaks = plot(tOG, Val_LNR, xlab = "t [sec]", ylab = "η [m]",  label = "Linear")
        plot!(L_PeakPos[1:NP], L_PeakVal[1:NP], seriestype=:scatter, ms=2, mc=:red, lab = "Peaks",legend=:bottomleft, legendcolumns=3)
        plot!([tOG[1];tOG[end]],[LMinPeakVal;LMinPeakVal], lab="$cstd η", mc=:red, ls=:dash)
        plot!(title = "Peaks of Linear signal", xlab = "t [sec]", ylab = "η [m]")

        # 8: Peak sections and mean of sections
        plta =  plot(tₑᵥ, events, legend = false, grid=true)
        plot!(tmean,NLN_mean,line = :dash, color = :red, linewidth = 2, label = "Mean")
        plot!(title = "Non-linear signal", xlab = "t [sec]", ylab = "η [m]")

        pltb = plot(Ltₑᵥ, Levents, legend = false, grid=true)
        plot!(tmeanLNR,LNR_mean,line = :dash, color = :red, linewidth = 2, label = "Mean")
        plot!(title = "Linear signal", xlab = "t [sec]", ylab = "η [m]")

        plt_ab = plot(plta, pltb, layout = @layout [a ; b ])

        # 9: Compare event spectrum with full spectrum
        fr_tmp, mag_tmp, _ = one_side_asp(Aₑₑ,tₑₑ)
        plt_evsp = plot(FR1[flb:fub], Â₀₀[flb:fub], line =:dashdot, lab="Non-linear")
        plot!(fr_tmp, mag_tmp, line =:solid, lab="Event")
        plot!(fⱼ, η̂, lab="JONSWAP", color=:red, lw=2)
        plot!(yscale=:log10, legendcolumns=1)
        plot!(xlim=(0,fcut), ylim=(1e-6,1), minorgrid=true)
        plot!(title = "Comparison of Spectra", xlab = L"f~[Hz]", ylab = L"S(f)~[m^2 s]")

        # 10: Illustrate position of events in the full signal
        h = PeakVal[evID]
        plt_groups = plot(tOG, A̅00, xlab = "t [sec]", ylab = "η [m]",  label = "Surface elevation")
        for i ∈ 1:min(20, length(PeakVal))
            local c = PeakPos[i]
            local w = 25
            plot!(Shape([c-w,c-w,c+w,c+w],[-h,h,h,-h]), color=:red, opacity=0.2)
        end
        plot!(legend = false)
        
        # 11: Illustrate position of selected event in the full signal
        c = PeakPos[evID]; w = 50
        plt_maxev = plot(tOG, A̅00, xlab = "t [sec]", ylab = "η [m]",  label = "Surface elevation")
        plot!(Shape([c-w,c-w,c+w,c+w],[-h,h,h,-h]), color=:red, opacity=0.3, label = "Event")
        plot!(xlim=(max(tOG[1],c-10*w), min(tOG[end],c+10*w)))

        # 12: Instantaneous frequency of event and event limits
        plt_ominst = plot(xlab="t [s]", ylab="ω [rad/s]", title="Event selection")
        plot!(tOG[lb-rng:ub+rng], omH[lb-rng:ub+rng], lab=L"ω_{inst}")
        plot!([tOG[PeakId[evID]];tOG[PeakId[evID]]],[minimum(omH[lb-rng:ub+rng]);maximum(omH[lb-rng:ub+rng])], lab="peak",ls=:dash)
        plot!([tOG[lb];tOG[lb]],[minimum(omH[lb-rng:ub+rng]);maximum(omH[lb-rng:ub+rng])], lab="lb",ls=:dash)
        plot!([tOG[ub];tOG[ub]],[minimum(omH[lb-rng:ub+rng]);maximum(omH[lb-rng:ub+rng])], lab="ub", ls=:dash)
        
        # 13: Interpolated event against original event
        plt_itp = plot(tEV[1:Lₚ], real(EV[1:Lₚ]), line = :solid, lw = 1, legend=:topleft, lab = "Linear (interp.)")
        plot!(twiny(), tpart, A̅00[lb:ub], line =:dash, color =:red, lw = 1, ylab = "η [m]", legend=:topright, lab = "Non-linear")
        plot!(title = "Event at Highest Peak: Original vs Interpolated", xlab = "t [sec]")

        # 14: Selected event & interpolation - Spectral analysis
        plt_mag_itp = plot(fr_ev, mag_ev, lw=2, ylab = L"Amplitude", lab="Event")
        plot!(FR, MAG, line =:dot, lw=2, lab="Interpolation")
        plot!(xlim=(0,fcut))

        plt_phi_itp = plot(xlab = "f [Hz]", ylab = L"\phi [rad]")
        plot!(fr_ev, unwrap(phi), lw=2, lab="Event")
        plot!(FR, PHI, lw=2, lab="Interpolation")
        plot!(xlim=(0,fcut))
        plot!(ylim=(minimum(PHI[1:icut]),maximum(PHI[1:icut])))

        plt_spitp = plot(plt_mag_itp, plt_phi_itp, layout = @layout [a; b])
    end

    # Plots to display
    display(pltϕ)
    display(plt_sp)
    display(plt1)
    display(plt_recon)
    display(plt_tcomp)
    display(plt_spcomp)
    
    if fevid
        display(plt_peaks)
        display(plt_Lpeaks)
        display(plt_ab)
        display(plt_evsp)
        display(plt_groups)
        display(plt_maxev)
        display(plt_ominst)
        display(plt_itp)
        display(plt_spitp)
    end

    # Figures to save
    fgnm = joinpath(DecFigs,"decomp_spectra")
    savefig(plt_sp, fgnm*".png"); savefig(plt_sp, fgnm*".svg")
    fgnm = joinpath(evdir,"decomp_elev")
    savefig(plt1, fgnm*".png");   savefig(plt1, fgnm*".svg")
    fgnm = joinpath(DecFigs,"spec_comp")
    savefig(plt_spcomp, fgnm*".png"); savefig(plt_spcomp, fgnm*".svg")

    if fevid
        fgnm = joinpath(DecFigs,"elev_peaks_nln")
        savefig(plt_peaks, fgnm*".png"); savefig(plt_peaks, fgnm*".svg")
        fgnm = joinpath(DecFigs,"elev_peaks_lin")
        savefig(plt_Lpeaks, fgnm*".png"); savefig(plt_Lpeaks, fgnm*".svg")
        fgnm = joinpath(DecFigs,"elev_pulses_compare")
        savefig(plt_ab, fgnm*".png"); savefig(plt_ab, fgnm*".svg")
        fgnm = joinpath(DecFigs,"events")
        savefig(plt_groups, fgnm*".png");   savefig(plt_groups, fgnm*".svg")
        fgnm = joinpath(evdir,"event_$evID")
        savefig(plt_maxev, fgnm*".png");    savefig(plt_maxev, fgnm*".svg")
        fgnm = joinpath(evdir,"sp_ev_$evID")
        savefig(plt_spitp, fgnm*".png");    savefig(plt_spitp, fgnm*".svg")
    end
end

#############################################################################################
if frec
    # Surface η (temporal)
    fid = joinpath(Decpath,"eta_lin") # File name
    open(fid, "w")
    writedlm(fid, [tOG Val_LNR], '\t')

    # Full surface elevation spectrum
    fid = joinpath(Decpath,"eta_spec") # File name
    open(fid, "w")
    writedlm(fid, [FR1 Â₀₀ ϕ₀₀], '\t')
    
    # 1st order elevation spectrum
    fid = joinpath(Decpath,"eta_lin_spec") # File name
    open(fid, "w")
    writedlm(fid, [FR1 AMPS[:,1] ARG[:,1]], '\t')

    # JONSWAP parameters
    fid = joinpath(Decpath,"JONSWAP_pars") # File name
    open(fid, "w")
    head = ["Hₛ" "Tₚ" "γ" "Tᵢ" "Tₑ"]
    row = round.([Hₛ Tₚ γ 1/fcut Tₑ]*1e3)./1e3
    writedlm(fid, [head; row], '\t')

    if fevid
        # Event
        fid = joinpath(evdir,"event") # File name
        open(fid, "w")
        writedlm(fid, [tₑₑ Aₑₑ], '\t')

        fid = joinpath(evdir,"event_lin") # File name
        open(fid, "w")
        writedlm(fid, [tₑₑ Aₑₑ¹], '\t')

        fid = joinpath(evdir,"SP_event") # File name
        open(fid, "w")
        writedlm(fid, [fr_ev mag_ev phi], '\t')

        fid = joinpath(evdir,"event_int") # File name
        open(fid, "w")
        writedlm(fid, [real(EV) imag(EV) angle.(EV)], '\t')

        fid = joinpath(evdir,"SP_event_int") # File name
        open(fid, "w")
        writedlm(fid, [FR MAG PHI], '\t')
    end
end
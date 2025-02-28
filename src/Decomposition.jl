module Decomposition
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

# Instantaneous frequency, Event Identification & Statistics
ωᴴ, evsep, ω̅ᵢₙₛₜ, ΔTᵉ, tΔT, idᵤₚ, Nᵤₚ, evAmax = instant_freq(A̅00,tOG)

#############################################################################################
# Find peaks above threshold
MinPeakDist = WAVEGEM.Tᵢ

## Non-linear surface eleveation
MinPeakVal = 2*std(A̅00) # Amplitude constraint (Hₛ/2)
PeakVal, PeakPos, PeakId = peaks_max_ext(A̅00, tOG, MinPeakVal, MinPeakDist,true)
Nqp = length(PeakVal)
println("Number of qualified peaks = $Nqp") 
Nqe = length(findall(evAmax .> MinPeakVal)) 
println("Number of qualified events = $Nqe")

## Linear surface eleveation
MinPeakVal = 2*std(Val_LNR) # Amplitude constraint (Hₛ/2)
L_PeakVal, L_PeakPos, L_PeakId = peaks_max_ext(Val_LNR, tOG, MinPeakVal, MinPeakDist,true)

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
    dfⱼ = abs(fⱼ[end]-fⱼ[1])/(length(fⱼ)-1)
    η̂ = sqrt.(2*Sⱼ*dfⱼ)

    # 1: Peaks on top of original signal
    plt_peaks = plot(tOG, A̅00, label = "Non-linear", palette=[cb[8];cb[11]])
    plot!(PeakPos, PeakVal, seriestype=:scatter, ms=2, lab = "Peaks",legend=:bottomleft, legendcolumns=2)
    plot!(title = "Peaks of Non-linear signal", xlab = "t [sec]", ylab = "η [m]")

    # 2: Peaks on top of linear signal
    plt_Lpeaks = plot(tOG, Val_LNR, xlab = "t [sec]", ylab = "η [m]",  label = "Linear", palette=[cb[8];cb[11]])
    plot!(L_PeakPos, L_PeakVal, seriestype=:scatter, ms=2, lab = "Peaks",legend=:bottomleft, legendcolumns=2)
    plot!(title = "Peaks of Linear signal", xlab = "t [sec]", ylab = "η [m]")

    # 3: Spectral decomposition
    flb = 1; fub = Int(round(fcut/df))

    plt_sp = plot(FR1[flb:fub], Â₀₀[flb:fub], line =:dashdot, lab="Non-linear", palette=[cb[8];cb[11];cb[5];cb[2];cb[9]])
    plot!(FR1[flb:fub], AMPS[flb:fub,1], line =:solid, lw=0.5, label = L"\mathbf{1_{st}}")
    plot!(FR1[flb:fub], AMPS[flb:fub,2], line =:solid, lw=0.5, label = L"\mathbf{2_{nd}^{-}}")
    plot!(FR1[flb:fub], AMPS[flb:fub,3], line =:dash, lw=0.5, label = L"\mathbf{2_{nd}^{+}}")
    plot!(FR1[flb:fub], AMPS[flb:fub,4], line =:dot, lw=0.5, label = L"\mathbf{3_{rd}}")
    plot!(fⱼ, η̂, lab="JONSWAP",  color=:red, lw=2)
    plot!(title = "Spectral Decomposition", xlab = L"f~[Hz]", ylab = L"S(f)~[m^2 s]")
    plot!(yscale=:log10, legendcolumns=1)
    plot!(xlim=(0,fcut), ylim=(1e-6,1), minorgrid=true)

    # 4: Compare spectra of input signal and OW3D simulation
    elev_fid::String = joinpath(rundir,"0","eta_t")
    open(elev_fid, "r")   
    elev = readdlm(elev_fid, '\t', Float64, skipstart=1, '\n')
    tinp, ηinp = elev[:,1], elev[:,2]
    frinp, maginp, _ = one_side_asp(ηinp,tinp)

    plt_comp = plot(FR1[flb:fub], Â₀₀[flb:fub], line =:dot, lab="Non-linear", palette=[cb[8];cb[11];cb[5]])
    plot!(FR1[flb:fub], AMPS[flb:fub,1], line =:dashdot, lw=1, opacity=0.8, lab = L"\mathbf{1_{st}}")
    plot!(frinp, maginp, lab="Linear Sea State")
    plot!(fⱼ, η̂, lab="JONSWAP", color=:red, lw=2)
    plot!(yscale=:log10, legendcolumns=1)
    plot!(xlim=(0,fcut), ylim=(1e-6,1), minorgrid=true)
    plot!(title = "Comparison of Spectra", xlab = L"f~[Hz]", ylab = L"S(f)~[m^2 s]")

    # 5: Illustrate events based on instantaneous frequency
    plt_evid = plot(xlab="t [s]", ylab="η [m]", title="Event identification based on instantaneous frequency", palette=[cb[8];cb[5]])
    plot!(tOG, A̅00, label = "Surface elevation")
    plot!(tOG, evsep, lab=:false)
    plot!(tOG[idᵤₚ], zeros(Float64, Nᵤₚ), seriestype=:scatter, ms=2, mc=:red, lab = "Event limits")

    # 6: Event durations
    plt_evdur = plot(bar(tΔT, ΔTᵉ, leg=:false, bar_width=ΔTᵉ, xlab="t [s]", ylab="ΔT [s]"), yaxis=:right) 
    plot!(twinx(), PeakPos, PeakVal, ylab="A [m]", seriestype=:scatter, ms=2, mc=:red, ylim=(minimum(PeakVal)-2,maximum(PeakVal)+2))
    title!("Event duration")    
    xlims!(0, tOG[end])

    # 7: Peak of each event
    plt_ADT = plot(xlab="ΔT [s]", ylab="A [m]", title="Event peak and duration")
    plot!(ΔTᵉ, evAmax, seriestype=:scatter, ms=2, mc=:red)
    plot!(range(0,maximum(ΔTᵉ),Nᵤₚ), 3*std(evAmax)*tanh.(1/(sqrt(2)/2*mean(ΔTᵉ)) * (range(0,maximum(ΔTᵉ),Nᵤₚ))))
    plot!(range(0,maximum(ΔTᵉ),Nᵤₚ), std(evAmax) .+ 3*std(evAmax)*tanh.(1/(sqrt(2)/2*mean(ΔTᵉ)) * (range(0,maximum(ΔTᵉ),Nᵤₚ))), line=:dot)
    plot!(range(0,maximum(ΔTᵉ),Nᵤₚ), 3*std(evAmax)*tanh.(1/(sqrt(2)/2*mean(ΔTᵉ)) * (range(0,maximum(ΔTᵉ),Nᵤₚ))) .- std(evAmax), line=:dot)

    mean(evAmax)
    # 8: Distribution of event durations (probably Rice)
    # sDT = std(ΔTᵉ) / sqrt(2-π/2)
    sDT = mean(ΔTᵉ) / sqrt(π/2)
    xDT = range(0,maximum(ΔTᵉ),Nᵤₚ)
    RayDT = xDT .* exp.((-xDT.^2)./(2*sDT^2)) ./ sDT^2
    # plot_DTdistr = histogram(ΔTᵉ, bins=Nᵤₚ) 
    plot_DTdistr = histogram(ΔTᵉ, normalize=:pdf) 
    plot!(xlab="ΔT [s]", ylab="P(ΔT)", title="Distribution of event durations (Rayleigh?)")
    plot!(xDT, RayDT)
    
    # Plots to display
    display(plt_peaks)
    display(plt_Lpeaks)
    display(plt_sp)
    display(plt_comp)
    display(plt_evid)
    display(plt_evdur)
    display(plt_ADT)
    display(plot_DTdistr)

    # Figures to save
    iname = ("elev_peaks_nln", "elev_peaks_lin", "decomp_spectra", "spect_comp")
    pname = (plt_peaks, plt_Lpeaks, plt_sp, plt_comp)

    for i ∈ eachindex(iname)
        ipath = joinpath(DecFigs,iname[i])
        savefig(pname[i], ipath*".png")
        savefig(pname[i], ipath*".svg")
    end
end

#############################################################################################
if frec
    # Surface η (temporal) - Nonlinear truncated signal
    fid_elev::String = "eta_nln" # File name
    open(joinpath(Decpath,fid_elev), "w")
    writedlm(joinpath(Decpath,fid_elev), [tOG A̅00], '\t')

    # Surface η (temporal) - 1st order component
    fid_elev::String = "eta_lin" # File name
    open(joinpath(Decpath,fid_elev), "w")
    writedlm(joinpath(Decpath,fid_elev), [tOG Val_LNR], '\t')

    # Surface η spectrum
    fid_elev::String = "eta_lin_spec" # File name
    open(joinpath(Decpath,fid_elev), "w")
    writedlm(joinpath(Decpath,fid_elev), [FR1 AMPS[:,1]], '\t')

    # Simulation parameters
    η̂ₘₐₓ = PeakVal[1]
    fid_elev::String = "sim_pars" # File name
    open(joinpath(Decpath,fid_elev), "w")
    head = ["Hₛ" "Tₚ" "γ" "Tᵢ" "Tₑ" "ηₘₐₓ" "T₂₋" "T₂₊" "ωᵢₙₛₜ" "Nqp" "Nqe"]
    row = [round.([Hₛ Tₚ γ WAVEGEM.Tᵢ Tₑ η̂ₘₐₓ T₂₋ T₂₊ ω̅ᵢₙₛₜ]*1e3)./1e3 [Int(Nqp) Int(Nqe)]]
    writedlm(joinpath(Decpath,fid_elev), [head; row], '\t')

    # OpenFAST elevation input file
    fid_elev::String = "Full.Elev" # File name
    open(joinpath(OFASTpath,"ExtElev","$(case_id)","$(run_id)",fid_elev), "w")
    writedlm(joinpath(OFASTpath,"ExtElev","$(case_id)","$(run_id)",fid_elev), [tOG Val_LNR], '\t')
end

if !frec && !fplot
    println("WARNING: Output recording and plots of Decomposition module have been suppressed (see Events module).")
end

end
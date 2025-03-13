module Events
using Plots, LaTeXStrings, DelimitedFiles, Dates, LinearAlgebra
using FFTW, DSP, Statistics, BSplineKit
import ColorSchemes.darkrainbow

gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))
cb = darkrainbow

# Load modules
include("WAVEGEM.jl")
import .WAVEGEM
include("Decomposition.jl")
import .Decomposition

# Include necessary scripts for functions
include("func/signal_processing.jl")
include("func/jonswap.jl")
include("func/peak_detect.jl")
include("func/directories.jl")
include("func/wave_theory.jl")
include("func/event_proc.jl")
include("func/text_process.jl")

#############################################################################################
# Import Global Variables & Module specific inputs
ρ, g, Ldom, d, _, _, γ, fcut = WAVEGEM.GlobInp0
pdir, run_id, phi_id, prb_id = WAVEGEM.GlobInp1
case_id, _, rundir = WAVEGEM.GlobPaths[2:4]
Decpath, _, DecFigs, postOFpath, _, OFASTpath = WAVEGEM.GlobPaths[9:end]
frec, fplot = WAVEGEM.Evflags
CET, evID = WAVEGEM.CET, WAVEGEM.evID
Wave = WAVEGEM.Wave

A̅00, Val_LNR = Decomposition.A̅00, Decomposition.Val_LNR
PeakId = Decomposition.PeakId
tOG, nₜ, dt = Decomposition.tOG, Decomposition.nₜ, Decomposition.dt
sig_comps, sig_recon = Decomposition.sig_comps, Decomposition.sig_recon
ωᴴ, idᵤₚ = Decomposition.ωᴴ, Decomposition.idᵤₚ
T₂₋, T₂₊ = Decomposition.T₂₋, Decomposition.T₂₊

#############################################################################################
# Make directories
if CET == 1        # Fairlead tension critical event
    case_str = "MaxFair"
elseif CET == 2    # Pitch critical event
    case_str = "MaxPitch"
elseif CET == 3    # CoM extreme displacement event 
    case_str = "MaxCoM"
else                # Response to extreme wave event
    case_str = "MaxWave"
end

if CET ∈ [1;2;3]
    f_tinst = joinpath(postOFpath,case_str*"_tinsts")   # OpenFAST events timestamps
    cont = parse_fxw(f_tinst, 0)                        # Read timestamps file
    tinst = cont[:,1]                                   # Vector of time instances
end

evdir = joinpath(Decpath,case_str,"EV$evID") # Output directory
make_dirs(2, joinpath(Decpath,case_str), evdir)

#############################################################################################
# Event Isolation & Analysis
## Event limits
if CET ∈ [1;2;3]
    lb = Int(tinst[evID]/dt+1)
else
    lb = PeakId[evID]
end
ub = lb+1

while sign(ωᴴ[lb]) == sign(ωᴴ[lb-1]) && lb > 2
    global lb = lb - 1
end

while sign(ωᴴ[ub]) == sign(ωᴴ[ub+1]) && ub < (nₜ-1)
    global ub = ub + 1
end
# Event exact duration
ΔTₑᵥ = tOG[ub] - tOG[lb]

# Extend event by one up- and one down-crossing prior to lb and after ub
while sign(Val_LNR[lb-1])-sign(Val_LNR[lb]) ≠ -2 && lb > 2
    global lb = lb - 1
end
while sign(Val_LNR[lb-1])-sign(Val_LNR[lb]) ≠ 2 && lb > 2
    global lb = lb - 1
end

while sign(Val_LNR[ub+1])-sign(Val_LNR[ub]) ≠ 2 && ub < (nₜ-1)
    global ub = ub + 1
end
while sign(Val_LNR[ub+1])-sign(Val_LNR[ub]) ≠ -2 && ub < (nₜ-1)
    global ub = ub + 1
end

tpart = tOG[lb:ub]; Lₚ = length(tpart)

Aₑₑ = A̅00[lb:ub]
Aₑₑ¹ = Val_LNR[lb:ub]
tₑₑ = tOG[lb:ub] .- tOG[lb]
Aᵐᵃˣ, idᵐᵃˣ = findmax(Aₑₑ) 
tᵐᵃˣ = tₑₑ[idᵐᵃˣ] .+ tOG[lb]

fr_ev, mag_ev, phi,_,_,Nfft, H = one_side_asp(Aₑₑ¹, tₑₑ)
itp_sp = interpolate(fr_ev, mag_ev, BSplineOrder(4))
itp_phi = interpolate(fr_ev, phi, BSplineOrder(4))
FR = range(0,fcut, Nfft)
MAG = itp_sp.(FR)
PHI = itp_phi.(FR)
EV = ifft(H)
tEV = zeros(Float64,Nfft)
tEV[1:Lₚ] = tₑₑ[:]

#############################################################################################
## PLOTS 
if fplot
    plt_mag_itp = plot(fr_ev, mag_ev, lw=2, ylab = L"Amplitude", lab="Event", palette=[cb[8];cb[11]])
    plot!(FR, MAG, line =:dot, lw=2, lab="Interpolation")
    plot!(xlim=(0,fcut))

    plt_phi_itp = plot(fr_ev, phi, lw=2, xlab = "f [Hz]", ylab = L"\phi [rad]", lab="Event", palette=[cb[8];cb[11]])
    plot!(FR, PHI, line =:dot, lw=2, lab="Interpolation")
    plot!(xlim=(0,fcut))

    plt_sp_itp = plot(plt_mag_itp, plt_phi_itp, layout = @layout [a; b])

    # 1: Decomposed signal and its components
    plt1 = plot(tpart, A̅00[lb:ub], lw = 2, lab = "Non-linear", palette=[cb[8];cb[11];cb[5];cb[2];cb[9]])
    plot!(tpart, Val_LNR[lb:ub], lw = 2, label = L"\mathbf{1_{st}}")
    plot!(tpart, sig_comps[lb:ub,2], line =:dot, lw = 2, lab = L"\mathbf{2_{nd}^{-}}")
    plot!(tpart, sig_comps[lb:ub,3], line =:dot, lw = 2, lab = L"\mathbf{2_{nd}^{+}}")
    plot!(tpart, sig_comps[lb:ub,4], line =:dot, lw = 2, lab = L"\mathbf{3_{rd}}")
    plot!(title = "Decomposed signal components", xlab = "t [sec]", ylab = "η [m]")
    # plot!(legendcolumns=3)

    # 2: Reconstructed signal against the original
    plt_recon = plot(tpart, A̅00[lb:ub], lw=2, label = "Non-linear",palette=[cb[8];cb[11]])
    plot!(tpart, sig_recon[lb:ub], line =:dash, lw=2, label = L"\sum components")
    plot!(title = "Original vs Reconstructed Signal", xlab = "t [sec]", ylab = "η [m]")

    # 3: Interpolated event against original event
    plt_itp = plot(tEV[1:Lₚ], real(EV[1:Lₚ]), lw = 2, ylab = "η [m]", lab = "Interpolation",  palette=[cb[8]], legend=:topright)
    plot!(twiny(), tpart, Aₑₑ¹, line =:dash, lw = 2, lab = L"1_{st}", palette=[cb[11]], legend=:bottomright)
    plot!(title = "Event - Linear SE", xlab = "t [sec]") 
    
    # 4: Highlight event position in full signal
    w = 5*(ub-lb)*dt
    h = maximum(A̅00[lb:ub])
    plt_maxev = plot(tOG, A̅00, xlab = "t [sec]", ylab = "η [m]",  label =:false, palette=[cb[8];cb[11]])
    plot!(Shape([tOG[lb],tOG[lb],tOG[ub],tOG[ub]],[-h,h,h,-h]), color=:red, opacity=0.3, label = "Event")
    plot!(xlim=(tOG[lb]-w, tOG[ub]+w))

    # 5: Selected event and event limits based on instantaneous frequency
    rng = Int(round(10/dt))
    if lb-rng < 1 || ub+rng > nₜ
        rng=Int(0)
    end
    plt_ominst = plot(xlab="t [s]", ylab="ω [rad/s]", title="Event instantaneous frequency", palette=[cb[11];cb[8];cb[4]])
    plot!(tOG[lb-rng:ub+rng], ωᴴ[lb-rng:ub+rng], lab=L"ω_{inst}", lw=2)
    plot!([tᵐᵃˣ; tᵐᵃˣ],[minimum(ωᴴ[lb-rng:ub+rng]);maximum(ωᴴ[lb-rng:ub+rng])], lab="peak",ls=:dash, lw=2)
    plot!([tOG[lb];tOG[lb]],[minimum(ωᴴ[lb-rng:ub+rng]);maximum(ωᴴ[lb-rng:ub+rng])], lab="lb",ls=:dash, lw=2)
    plot!([tOG[ub];tOG[ub]],[minimum(ωᴴ[lb-rng:ub+rng]);maximum(ωᴴ[lb-rng:ub+rng])], lab="ub", ls=:dash, lw=2)

    # 6: 2nd- component spectrum
    plt_2ndm_a = plot(xlab = "t [sec]", ylab = "η [m]", palette=[cb[5]], title=L"2_{nd}^{-}"*" component")
    plot!(tpart, sig_comps[lb:ub,2], line =:dot, lw = 2, lab=:false)

    _, mag_2nd,phi_2nd,_ = one_side_asp(sig_comps[lb:ub,2], tₑₑ)
    plt_2ndm_b = plot(xlab = "f [Hz]", ylab = L"S(f)~[m]", lw=2, xlim=(0,fcut), palette=[cb[8];cb[11];cb[5];cb[2]])
    plot!(fr_ev[2:end], mag_2nd[2:end], lw=2, lab=:false, xscale=:log10, xlim=(fr_ev[2],1))
    plot!([1/T₂₋; 1/T₂₋], [1e-12; maximum(mag_2nd)], line=:dot, lw=1, lab=L"f_2^{-}")
    plot!([1/T₂₊; 1/T₂₊], [1e-12; maximum(mag_2nd)], line=:dot, lw=1, lab=L"f_2^{-}")

    plt_2ndm_c = plot(xlab = "f [Hz]", ylab = L"\phi~[rad]", palette=[cb[8]])
    plot!(fr_ev[2:end], phi_2nd[2:end], lw=2, lab=:false, xscale=:log10)
    plot!(xlim=(fr_ev[2],1))

    plt_2ndm = plot(plt_2ndm_a, plt_2ndm_b, plt_2ndm_c, layout = @layout [a; b c])
    
    display(plt1)
    display(plt_recon)
    display(plt_itp)
    display(plt_maxev)
    display(plt_ominst)
    display(plt_sp_itp)
    display(plt_2ndm)

    # Figures to save
    iname = ("decomp_elev", "elev_compare", "XEvent", "SP_event", "om_inst", "sp_2nd-")
    pname = (plt1, plt_recon, plt_maxev, plt_sp_itp, plt_ominst, plt_2ndm)

    for i ∈ eachindex(iname)
        ipath = joinpath(evdir,iname[i])
        savefig(pname[i], ipath*".png")
        savefig(pname[i], ipath*".svg")
    end
end

#############################################################################################
if frec
    # Event
    fid_ev = joinpath(evdir,"event") # File name
    open(fid_ev, "w")
    writedlm(fid_ev, [tₑₑ Aₑₑ], '\t')

    fid_ev = joinpath(evdir,"event_lin") # File name
    open(fid_ev, "w")
    writedlm(fid_ev, [tₑₑ Aₑₑ¹], '\t')

    # Event 2nd- component (temporal)
    fid_ev = joinpath(evdir,"eta_2nd-") # File name
    open(fid_ev, "w")
    writedlm(fid_ev, [tₑₑ sig_comps[lb:ub,2]], '\t')

    fid_ev = joinpath(evdir,"SP_event") # File name
    open(fid_ev, "w")
    writedlm(fid_ev, [fr_ev mag_ev phi], '\t')

    # 2nd- component spectrum
    fid_ev = joinpath(evdir,"eta_2nd-_spec") # File name
    open(fid_ev, "w")
    writedlm(fid_ev, [fr_ev mag_2nd], '\t')

    fid_ev = joinpath(evdir,"event_int") # File name
    open(fid_ev, "w")
    writedlm(fid_ev, [real(EV) imag(EV) angle.(EV)], '\t')

    fid_ev = joinpath(evdir,"SP_event_int") # File name
    open(fid_ev, "w")
    writedlm(fid_ev, [FR MAG PHI], '\t')

    # Event parameters
    Hᵉₛ = 4*std(A̅00[lb:ub])
    η̂ᵉₘₐₓ = maximum(A̅00[lb:ub])
    Tᵉₚ =  1/FR[findmax(MAG)[2]]
    Tᵉ₂₋ = 1/fr_ev[findmax(mag_2nd)[2]]

    fid_ev = joinpath(evdir,"ev_pars") # File name
    open(fid_ev, "w")
    head = ["ηₘₐₓ" "Hₛ" "Tₚ" "T₂₋" "ΔTₑᵥ"]
    row = round.([η̂ᵉₘₐₓ Hᵉₛ Tᵉₚ Tᵉ₂₋ ΔTₑᵥ]*1e3)./1e3
    writedlm(fid_ev, [head; row], '\t')
end
    
end
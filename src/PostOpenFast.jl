module PostOpenFast
using Plots, LaTeXStrings, DelimitedFiles, Dates, LinearAlgebra
using FFTW, DSP, Statistics, BSplineKit

include("WAVEGEM.jl")
import .WAVEGEM

# Include necessary scripts for functions
include("func/signal_processing.jl")
include("func/text_process.jl")
include("func/peak_detect.jl")
include("func/jonswap.jl")
include("func/event_proc.jl")

import ColorSchemes.darkrainbow
gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))
cb = darkrainbow
#############################################################################################
# Import Global Variables & Module specific inputs
ρ, g, _, d, _, _, γ, fcut = WAVEGEM.GlobInp0
runstr, case_id, casedir, rundir, phipath, OW3Dcdir, OW3Drdir, OW3Dphipath, 
            Decpath, DecEvs, DecFigs, postOFpath, postOW3Dpath, OFASTpath = WAVEGEM.GlobPaths
Wave = WAVEGEM.Wave
full, fev = WAVEGEM.POFflags
CET, evID = WAVEGEM.CET, WAVEGEM.evID
T₀ₛ, T₀ₕ, T₀ₚ = WAVEGEM.FOWT[:]

#############################################################################################
# Open and read files
## Full simulation
OFASTout = joinpath(postOFpath,"outD_")   # OpenFAST output file
cont = parse_fxw(OFASTout, 5)
heads = readdlm(joinpath(postOFpath,"dataHdr"))

## Surface elevations for comparative plots
sefname = joinpath(Decpath,"eta_nln")      # Surface elevation from OW3D (0 phase shift)
ηOW3D = parse_fxw(sefname, 0)

#############################################################################################
# Assign Variables
NoRows = size(cont)[2]

# Nₜ = size(cont)[1]
# dt = 0.01
# t = zeros(Float64,Nₜ)
# [t[i] = (i-1)*dt for i ∈ 1:Nₜ]  
# Surge, Heave, Pitch, Wave1Elev, Wave1Elev1, FAIRTEN1, FAIRTEN2, 
# HydroFxi, HydroFzi, HydroMyi, PtfmTAxt, PtfmTAzt = [cont[:,i] for i = 1:NoRows]

t, Surge, Heave, Pitch, Wave1Elev, Wave1Elev1, FAIRTEN1, FAIRTEN2, 
HydroFxi, HydroFzi, HydroMyi, PtfmTAxt, PtfmTAzt = [cont[:,i] for i = 1:NoRows]
aCoM = sqrt.(PtfmTAxt.^2 .+ PtfmTAzt.^2)  # CoM acceleration
Nₜ = length(t)       
dt = t[2]-t[1]

#############################################################################################
# Critical Event statistics & selection
RespVar = zeros(Float64,Nₜ) # The response variable to be analysed

if CET == 1        # Fairlead tension critical event
    RespVar[:] = FAIRTEN2[:]
    case_str = "MaxFair"
elseif CET == 2    # Pitch critical event
    RespVar[:] = Pitch[:]
    case_str = "MaxPitch"
elseif CET == 3    # CoM extreme displacement event 
    RespVar[:] = disp[:]
    case_str = "MaxCoM"
else                # Response to extreme wave event
    RespVar[:] = Wave1Elev[:]
    case_str = "MaxWave"
end
figOF = joinpath(postOFpath,case_str)
## Make directory for specified event if it doesn't exist
if !isdir(figOF)
    mkdir(figOF)
end

# Find all maxima above a given threshold of the selected response time history
PosPeaks, _ = peaks_max_ext(RespVar,t, RespVar[1], 0, false)

## Find and sort extreme maxima in the selected response time history
MinPeakVal = 5*std(RespVar) .+ mean(RespVar)    # Amplitude threshold (default=3*std)
Aᵢ, tᵢ, i⁺, _ = peaks_max_ext(RespVar,t, MinPeakVal, 0, true)
Nce = length(Aᵢ)        # Number of peaks above threshold

# Remove duplicate events
## Instantaneous frequency
ωᴴ,_ = instant_freq(RespVar,t)
## Start point of corresponding event
lbs = Array{Int64}(undef,0)
for i ∈ 1:Nce
    ilb = i⁺[i]
    while sign(ωᴴ[ilb]) == sign(ωᴴ[ilb-1]) && ilb > 2
        ilb = ilb - 1
    end
    push!(lbs,ilb)
end
duplicates, non_duplicates = separate_duplicates(lbs)
println("Duplicate events start time= ", t[Int.(duplicates)])

### Store time instances of maxima
f_tinst = joinpath(postOFpath,case_str*"_tinsts")
open(f_tinst, "w")
writedlm(f_tinst, [tᵢ Aᵢ], '\t')


# Time histories
if fev
    # Event
    cont_ev = parse_fxw(figOF*"/outD_EV$evID", 5)
    # FWG
    cont_FWG = parse_fxw(figOF*"/outD_EV$(evID)_ReFoGWs", 5)
    # # DAM
    # cont_DAM = parse_fxw_pf(figOF*"/outD_DAM", 0, 5)
    # # 2AM
    # cont_2AM = parse_fxw_pf(figOF*"/outD_2AM", 0, 5)
    # # ALT_2AM
    # cont_ALT_2AM = parse_fxw_pf(figOF*"/outD_ALT_2AM", 0, 5)
    # # SFWG
    # # cont_SFWG = parse_fxw_pf(figOF*"/outD_SFWG", 0, 5)
    # # DWG
    # cont_DWG = parse_fxw_pf(figOF*"/outD_DWG", 0, 5)
    # # BEAT
    # cont_BEAT = parse_fxw_pf(figOF*"/outD_BEAT", 0, 5)
    # # NewWave
    # cont_NW = parse_fxw_pf(postOFpath*"/outD_NW", 0, 5)
    trange = Int(round(cont_ev[end,1]/2/dt))
    imax = findmax(cont_ev[:,8])[2] # Position of max fairten  
else
    trange = Int(round(2^7/dt))
    imax = 0
end

# Truncated time vector for event plots 
peak = [tᵢ[evID];i⁺[evID]]    # Time and interval of event peak
lb = Int(peak[2]) - imax
ub = Int(lb+2*trange-1) 
tₑᵥ = t[lb:ub];     Lₑᵥ = length(tₑᵥ)

#############################################################################################
# PLOTS
## Plot all stored variables of the whole simulation
if full
    # Plot OpenFast results
    for i ∈ 2:NoRows
        plti = plot(xlab = "t [s]", title = "$(heads[i])", legend=:topleft, palette=[cb[11]])
        plot!(t, cont[:,i], lw=2, lab = "Full sim")
        display(plti)
        savefig(joinpath(postOFpath,"$(heads[i]).svg"))
        savefig(joinpath(postOFpath,"$(heads[i]).png"))
    end

    # Wave Spectra
    FR, MAG, _ = one_side_asp(Wave1Elev, t)               # OpenFAST - spectrum of surface elevation
    FR_1st, MAG_1st, _ = one_side_asp(Wave1Elev1, t)      # OpenFAST - spectrum of 1st order elevation
    FR_2nd, MAG_2nd, _ = one_side_asp(Wave1Elev .- Wave1Elev1, t)   # OpenFAST - spectrum 2nd order
    FR_OW3D, MAG_OW3D, _ = one_side_asp(ηOW3D[:,2], ηOW3D[:,1])     # Spectrum of OW3D elevation (nln)
    fₚ = FR[findmax(MAG)[2]]

    # Plot OpenFAST surface elevation spectra
    plt = plot(xlab = L"f~[Hz]", ylab = L"S(f)~[m]", palette=[cb[11];cb[4];cb[8]]) 
    plot!(FR, MAG, yscale=:log10, line =:dot, lab = "OpenFAST")
    plot!(FR_1st, MAG_1st, yscale=:log10, opacity=0.8, lab = L"1^{st}~\mathrm{order}")
    plot!(FR_2nd, MAG_2nd, yscale=:log10, opacity=0.6, lab = L"2^{nd}~\mathrm{order}")
    plot!(xlim=(0,fcut), ylim=(1e-6,1),  minorgrid=true)
    display(plt)
    savefig(joinpath(postOFpath,"OFWaveSpectra.svg"))
    savefig(joinpath(postOFpath,"OFWaveSpectra.png"))

    # Comparative plots of surface elevation spectra: OpenFAST vs OW3D 
    plt = plot(xlab = L"f~[Hz]", ylab = L"S(f)~[m]", legendcolumns=:2, palette=[cb[11];cb[4];cb[8]]) 
    plot!(FR_OW3D, MAG_OW3D, lab = "OW3D")
    plot!(FR, MAG, line =:dash, lw = 0.5, lab = "OpenFAST")
    plot!(1/T₀ₛ*ones(Float64,100), maximum(MAG)*range(0,1,100), line=:dash,  lw = 2, lab = L"f_{0}^{surge}")
    plot!(1/T₀ₕ*ones(Float64,100), maximum(MAG)*range(0,1,100), line=:dash, lw = 2, lab = L"f_{0}^{heave}")
    plot!(1/T₀ₚ*ones(Float64,100), maximum(MAG)*range(0,1,100), line=:dash, lw = 2, lab = L"f_{0}^{heave}")
    plot!(xlim=(0, 0.625))
    display(plt)
    savefig(joinpath(postOFpath,"CompWaveSpectra.svg"))
    savefig(joinpath(postOFpath,"CompWaveSpectra.png"))

    # Distribution of RespVal peaks
    plt_PP = histogram(PosPeaks, normalize=:pdf)
    display(plt_PP)

    # Plot the selected response variable
    pltRV = plot(t, RespVar, xlab="t [s]", ylab="RespVar", title=case_str, legend=:false)
    plot!([t[1],t[end]], [MinPeakVal,MinPeakVal], line=:dot)
    display(pltRV)
end

## Comparative multi-plot of key responses to selected extreme event
### Subplot 1 - Surge and CoM displacement
plt_eva = plot(legend=:topleft, palette=[cb[10];cb[10]])
plot!(tₑᵥ, Surge[lb:ub], lw=2, ylab="[m]", lab="Surge")
# plot!(tₑᵥ, aCoM[lb:ub], lw=1, lab=L"a_{CM}", line=:dot)
plot!(tₑᵥ, maximum(Surge)*ones(Float64,Lₑᵥ), lab="max", line=:dash)
# plot!(tₑᵥ, maximum(aCoM)*ones(Float64,Lₑᵥ), lab="max", line=:dot, lw=1)
### Subplot 2 - Pitch angle
plt_evb = plot(legend=:topleft, palette=[cb[4];cb[4]])
plot!(tₑᵥ, Pitch[lb:ub], lw=2, ylab="[deg]", lab="Pitch")
plot!(tₑᵥ, maximum(Pitch)*ones(Float64,Lₑᵥ), lab="max", line=:dash)
### Subplot 3 - Fairlead tension
plt_evc = plot(xlab="t [s]", legend=:topleft, palette=[cb[9];cb[9]])
plot!(tₑᵥ, FAIRTEN2[lb:ub], lw=2, ylab="[kN]", lab = "Fairlead")
plot!(tₑᵥ, maximum(FAIRTEN2)*ones(Float64,Lₑᵥ), lab="max", line=:dash)
### Subplot 4 - Surface elevation and Heave response
plt_evd = plot(xlab="t [s]", legend=:bottomleft, legendcolumns=:2, palette=[cb[1];cb[8];cb[1];cb[8]])
plot!(tₑᵥ, Wave1Elev[lb:ub], lw=2,  lab=L"\eta")
# plot!(tₑᵥ, Wave1Elev1[lb:ub], lw=2,  lab=L"\eta^1")
plot!(tₑᵥ, Heave[lb:ub], lw=2, ylab="[m]", lab="Heave")
plot!(tₑᵥ, maximum(Wave1Elev)*ones(Float64,Lₑᵥ), lab=L"max(\eta)", line=:dash)
# plot!(tₑᵥ, maximum(Wave1Elev1)*ones(Float64,Lₑᵥ), lab=L"max(\eta^1)", line=:dash)
plot!(tₑᵥ, maximum(Heave)*ones(Float64,Lₑᵥ), lab=L"max(z)", line=:dash)
### Concentrated plot (1+2+3+4)
plt_event = plot(plt_eva, plt_evb, plt_evc, plt_evd, layout = @layout [a b; c d])
display(plt_event)
savefig(joinpath(figOF,"EE$(evID)_resp.svg"))
savefig(joinpath(figOF,"EE$(evID)_resp.png"))


#############################################################################################
if fev
    for i ∈ 2:NoRows
        plti = plot(xlab = "t [s]", title = "$(heads[i])", legend=:topleft, palette=[cb[11];cb[4];cb[8];cb[1]])
        plot!(tₑᵥ, cont[lb:ub,i], lw=3, lab = "Full sim")
        plot!(tₑᵥ, cont_ev[:,i], lw=2, lab = "Event")
        plot!(tₑᵥ, cont_FWG[:,i], line=:dot, lw=2, lab = L"ReFoGWs")
        # plot!(tₑᵥ, cont_SFWG[:,i], line=:dot, lw=2, lab = L"\sum g_n (t)")
        # plot!(tₑᵥ, cont_NW[:,i], line=:dot, lw=2, lab = "Focused Wave")
        display(plti)
        # savefig(figOF*"/$(heads[i]).svg")
        # savefig(figOF*"/$(heads[i]).png")
    end

    pltFairA = plot(xlab = "t [s]", title = "$(heads[8])", legend=:topleft, palette=[cb[11];cb[4];cb[8];cb[1]])
    plot!(tₑᵥ, FAIRTEN2[lb:ub], lw=3, lab = "Full sim")
    plot!(tₑᵥ, cont_FWG[:,8], lw=2, lab = L"ReFoGWs")
    display(pltFairA)
    savefig(joinpath(figOF,"FairA_ev$evID.svg"));  savefig(joinpath(figOF,"FairA_ev$evID.png"))
end
end
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
CEid = WAVEGEM.CEid
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
Nₜ = length(t)       
dt = t[2]-t[1]


aCoM = sqrt.(PtfmTAxt.^2 .+ PtfmTAzt.^2)  # CoM acceleration

#############################################################################################
# Critical Event statistics & selection
RespVar = zeros(Float64,Nₜ) # The response variable to be analysed

if CEid == 1        # Fairlead tension critical event
    RespVar[:] = FAIRTEN2[:]
    case_str = "MaxFair"
elseif CEid == 2    # Pitch critical event
    RespVar[:] = Pitch[:]
    case_str = "MaxPitch"
elseif CEid == 3    # CoM extreme displacement event 
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

## Find and sort maxima in the selected response time history
MinPeakVal = 3*std(RespVar) .+ mean(RespVar)    # Amplitude threshold (default=3*std)
Aᵢ, tᵢ, i⁺, _ = peaks_max_ext(RespVar,t, MinPeakVal, 0, true)
Nce = length(Aᵢ)        # Number of peaks above threshold

### Store time instances of maxima
f_tinst = joinpath(postOFpath,case_str*"_tinsts")
if !isfile(f_tinst)
    open(f_tinst, "w")
    writedlm(f_tinst, [tᵢ Aᵢ], '\t')
end

peak = [tᵢ[2];i⁺[2]]    # Time and interval of event peak

# Truncated time vector for event plots
trange = Int(round(255.5/dt));  
lb = Int(peak[2]-trange);      ub = Int(lb+2*trange)
# lb = Int(round(lb/10)*10);  ub = Int(round(ub/10)*10) 
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
savefig(joinpath(figOF,"EE_resp.svg"))
savefig(joinpath(figOF,"EE_resp.png"))

end
# Load necessary modules
using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles
using FFTW, DSP, Statistics, BSplineKit
import ColorSchemes.darkrainbow

gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))
cb = darkrainbow

# Include necessary scripts for functions
include("signal_processing.jl")
include("text_process.jl")
include("peak_detect.jl")
include("jonswap.jl")
include("directories.jl")

#############################################################################################
## Read from text files in library folder
#############################################################################################
# Directories
pdir::String = "library/SE/"    # Parent directory
cdir::String = "LCH"            # Case directory
rdir::String = "/1"             # Run directory
phidir::String = "/0"           # Phase shift corresponding directory
OW3Dprb::String = "/x1000"    # PostOW3D probe corresponding directory
# Files
OFASTout::String = "/outD_"   # OpenFAST output file
evfile::String = "event_1"      # Event file name
#############################################################################################
# Path definition 
rpath = pdir*cdir*rdir  # Path to run directory
phipath = rpath*phidir  # Path to phase shift corresponding directory
Dpath = rpath*"/Decomposition"     # Path to Decomposition directory
evdpath = Dpath*"/events"          # Path to events directory
postOW3Dpath = phipath*"/postOW3D"*OW3Dprb    # Path to OW3D postprocessing directory
postOFpath = phipath*"/postOFAST"     # Path to OpenFAST postprocessing directory
OFoutpath = postOFpath*OFASTout       # Path to OpenFAST output file
evfpath = evdpath*evfile           # Path to event file
figOF = postOFpath*"/"          # Path to OpenFast event output folder

# Flags
full = Bool(0)  # Flag to plot and record full simulation
fev = Bool(1)   # Flag to plot and record event simulations
EEid = 0          # 0: Max Fairlead tension, 1: Max Pitch, 2: Max Wave, 3: Max CoM displacement

# Open and read files
## Full simulation
cont = parse_fxw_pf(OFoutpath, 0, 5)
heads = readdlm(postOFpath*"/dataHdr")

#############################################################################################
fcut = 0.625;   dt = 0.01
# Semi-submersible FOWT eigenfrequencies
T₀ₛ¹ = 250.0000;    T₀ₛ² = 35.714107  
T₀ₕ = 17.543772;    T₀ₚ = 27.026892

Nₜ = size(cont)[1];     NoRows = size(cont)[2]
t = zeros(Float64,Nₜ);      [t[i] = (i-1)*dt for i ∈ 1:Nₜ]

Surge = cont[:,1];      Heave = cont[:,2];      Pitch = cont[:,3]
Wave1Elev = cont[:,4];  Wave1Elev1 = cont[:,5]; FAIRTEN2 = cont[:,7]

if full
    for i ∈ 1:12
        plti = plot(xlab = "t [s]", title = "$(heads[i])", legend=:topleft, palette=[cb[11]])
        plot!(t, cont[:,i], lw=2, lab = "Full sim")
        display(plti)
        savefig(figOF*"$(heads[i]).svg")
        savefig(figOF*"$(heads[i]).png")
    end

    # Wave Spectra
    ηOW3D = parse_fxw_pf(postOW3Dpath*"/eta", 0, 0)
    ηlin = parse_fxw_pf(Dpath*"/eta_lin", 0, 0)
    ηOG = parse_fxw_pf(Dpath*"/eta0", 0, 0)

    FRof, Mof, _ = one_side_asp(Wave1Elev, t)
    FRlin, Mlin, _ = one_side_asp(ηlin[:,2], ηlin[:,1])
    FRog, Mog, _ = one_side_asp(ηlin[:,2], ηlin[:,1])
    FR, MAG, _ = one_side_asp(Wave1Elev .- Wave1Elev1, t)

    fₚ = FRof[findmax(Mof)[2]]

    plt = plot(FRof, Mof, yscale=:log10, line =:dot, xlab = L"f~[Hz]", ylab = L"S(f)~[m^2 s]", lab = "OpenFAST", palette=[cb[11];cb[4];cb[8]])
    plot!(FRlin, Mlin, yscale=:log10, opacity=0.8, lab = "Linear")
    plot!(FR, MAG, yscale=:log10, opacity=0.6, lab = "2nd order")
    plot!(xlim=(0,fcut), ylim=(1e-6,1),  minorgrid=true)
    display(plt)
    savefig(figOF*"/OFWaveSpectra.svg")
    savefig(figOF*"/OFWaveSpectra.png")

    plt = plot(FRog, Mog, xlab = L"f~[Hz]", ylab = L"S(f)~[m^2 s]", lab = "Original", legendcolumns=:2, palette=[cb[11];cb[4];cb[8]])
    plot!(FRof, Mof, line =:dash, lw = 0.5, lab = "OpenFAST")
    plot!(1/T₀ₛ¹*ones(Float64,100), maximum(Mof)*range(0,1,100), line=:dash,  lw = 2, lab = L"f_{01}^{surge}")
    plot!(1/T₀ₛ²*ones(Float64,100), maximum(Mof)*range(0,1,100), line=:dash, lw = 2, lab = L"f_{02}^{surge}")
    plot!(1/T₀ₕ*ones(Float64,100), maximum(Mof)*range(0,1,100), line=:dash, lw = 2, lab = L"f_{0}^{heave}")
    plot!(1/T₀ₚ*ones(Float64,100), maximum(Mof)*range(0,1,100), line=:dash, lw = 2, lab = L"f_{0}^{heave}")
    plot!(xlim=(0, 0.625))
    display(plt)
    savefig(figOF*"/CompWaveSpectra.svg")
    savefig(figOF*"/CompWaveSpectra.png")
end

#############################################################################################
# Extreme Event selection
r = 41.216
θ₀ = -(π/2+atan(40.868/5.3412))
u = r*(cos.(θ₀.+Pitch*π/180).-cos(θ₀))
w = r*(sin.(θ₀.+Pitch*π/180).-sin(θ₀))
Fdisp = sqrt.((Surge .+ u).^2 .+ (Heave .+ w).^2)
disp = sqrt.(Surge.^2 .+ Heave.^2)

if EEid == 0
    peak = findmax(FAIRTEN2)
    figOF = figOF*"MaxFair"
elseif EEid == 1
    peak = findmax(Pitch)
    figOF = figOF*"MaxPitch"
elseif EEid == 2
    peak = findmax(Wave1Elev)
    figOF = figOF*"MaxWave"
elseif EEid == 3 
    peak = findmax(disp)
    figOF = figOF*"MaxCoM"

    plt_disp = plot(xlab="t~[s]", ylab="Displacement [m]", palette=:darkrainbow)
    plot!(t, disp, lab="CoM")
    plot!(t, Fdisp, lab="Fairled")
    display(plt_disp)
end

if !isdir(figOF)
    mkdir(figOF)
end

# Truncated time vector for event plots
trange = Int64(round(50/dt));  
lb = peak[2]-3*trange;      ub = lb+4*trange 
lb = Int(round(lb/10)*10);  ub = Int(round(ub/10)*10) 
tₑᵥ = t[lb:ub];     Lₑᵥ = length(tₑᵥ)

# Plot full simulation response to selected extreme event
plt_eva = plot(legend=:topleft, palette=[cb[10];cb[1]])
plot!(tₑᵥ, Surge[lb:ub], lw=2, ylab="[m]", lab="Surge")
plot!(tₑᵥ, disp[lb:ub], lw=1, lab="CoG", line=:dot)
plot!(tₑᵥ, maximum(Surge)*ones(Float64,Lₑᵥ), lab="max(x)", line=:dash)
plot!(tₑᵥ, maximum(disp)*ones(Float64,Lₑᵥ), lab="max(d)", line=:dot, lw=1)

plt_evb = plot(legend=:topleft, palette=[cb[5];cb[1]])
plot!(tₑᵥ, Pitch[lb:ub], lw=2, ylab="[deg]", lab="Pitch")
plot!(tₑᵥ, maximum(Pitch)*ones(Float64,Lₑᵥ), lab="max", line=:dash)

plt_evc = plot(xlab="t [s]", legend=:topleft, palette=[cb[7];cb[1]])
plot!(tₑᵥ, FAIRTEN2[lb:ub], lw=2, ylab="[kN]", lab = "Fairlead")
plot!(tₑᵥ, maximum(FAIRTEN2)*ones(Float64,Lₑᵥ), lab="max", line=:dash)

plt_evd = plot(xlab="t [s]", legend=:bottomleft, legendcolumns=:3, palette=[cb[11];cb[4];cb[8]])
plot!(tₑᵥ, Wave1Elev[lb:ub], lw=2,  lab=L"\eta")
plot!(tₑᵥ, Wave1Elev1[lb:ub], lw=2,  lab=L"\eta^1")
plot!(tₑᵥ, Heave[lb:ub], lw=2, ylab="[m]", lab="Heave")
plot!(tₑᵥ, maximum(Wave1Elev)*ones(Float64,Lₑᵥ), lab=L"max(\eta)", line=:dash)
plot!(tₑᵥ, maximum(Wave1Elev1)*ones(Float64,Lₑᵥ), lab=L"max(\eta^1)", line=:dash)
plot!(tₑᵥ, maximum(Heave)*ones(Float64,Lₑᵥ), lab=L"max(z)", line=:dash)

plt_event = plot(plt_eva, plt_evb, plt_evc, plt_evd, layout = @layout [a b; c d])
display(plt_event)
savefig(figOF*"/EE_resp.svg")
savefig(figOF*"/EE_resp.png")

#############################################################################################
# Time histories
if fev
    # Event
    cont_ev = parse_fxw_pf(figOF*"/outD_event", 0, 5)
    # Gaussian interpolation
    cont_G = parse_fxw_pf(figOF*"/outD_etaG", 0, 5)
    # Sum of EWGs
    cont_EWGsum = parse_fxw_pf(figOF*"/outD_EWGsum", 0, 5)
    # Focused Wave
    cont_NW = parse_fxw_pf(figOF*"/outD_etaFCSD", 0, 5)

    for i ∈ 1:12
        plti = plot(xlab = "t [s]", title = "$(heads[i])", legend=:topleft, palette=[cb[11];cb[4];cb[8];cb[1]])
        plot!(tₑᵥ, cont[lb:ub,i], lw=3, lab = "Full sim")
        plot!(tₑᵥ, cont_ev[:,i], lw=2, lab = "Event")
        plot!(tₑᵥ, cont_G[:,i], line=:dot, lw=2, lab = L"G(t)")
        plot!(tₑᵥ, cont_EWGsum[:,i], line=:dot, lw=2, lab = L"\sum g_n (t)")
        # plot!(tₑᵥ, cont_NW[:,i], line=:dot, lw=2, lab = "Focused Wave")
        display(plti)
        savefig(figOF*"/$(heads[i]).svg")
        savefig(figOF*"/$(heads[i]).png")
    end
end

# # Plot RAOs
# FR = readdlm("library/SE/LAF/2/0/postOFAST/FREQ", skipstart=5, '\n');
# PSD = parse_fxw_pf("SE/LAF/2/0/postOFAST/PSDs", 0, 5);

# pPSDs = plot(FR, PSD[:,1], xlab = L"f~[Hz]", ylab = L"PSD~[Unit^2/Hz]", lab = L"Surge", lw=3, yscale=:log10)
# plot!(FR, PSD[:,2], lab = L"Heave", lw=3, yscale=:log10)
# plot!(FR, PSD[:,3],  lab = L"Pitch", lw=3, yscale=:log10)
# plot!([fₚ, fₚ], [1e-7, findmax(PSD[:,3])[1]], line=:dash, lw=2, color=:red, lab = L"f_p", yscale=:log10)
# plot!(xlim=(0, 0.2), ylim=(1e-7, 1e3), title = "Power Spectral Densities of main DOFs", minorgrid=true)
# display(pPSDs)
# savefig(postOFpath*"/PSDs.svg")
# savefig(postOFpath*"/PSDs.png")

# T = 1 ./FR[2:end]
# pPSDs = plot(T, PSD[2:end,1], xlab = L"T~[s]", ylab = L"PSD~[Unit^2/Hz]", lab = L"Surge", lw=3, yscale=:log10)
# plot!(T, PSD[2:end,2], lab = L"Heave", lw=3, yscale=:log10)
# plot!(T, PSD[2:end,3],  lab = L"Pitch", lw=3, yscale=:log10)
# plot!([1/fₚ, 1/fₚ], [1e-7, findmax(PSD[:,3])[1]], line=:dash, lw=2, color=:red, lab = L"f_p", yscale=:log10)
# plot!(xlim=(0, 80), ylim=(1e-7, 1e3), title = "Power Spectral Densities of main DOFs", minorgrid=true)
# display(pPSDs)
# savefig(postOFpath*"/PSDs_T.svg")
# savefig(postOFpath*"/PSDs_T.png")


# Wave1Elev = cont[:,4]
# Wave1Elev2 = cont[:,5]
# Wave1Elev1 = cont[:,6]
# HydroFxi = cont[:,7]
# HydroFyi = cont[:,8] 
# HydroFzi = cont[:,9] 
# HydroMxi = cont[:,10] 
# HydroMyi = cont[:,11]  
# HydroMzi = cont[:,12]
# FAIRTEN1 = cont[:,13] 
# FAIRTEN2 = cont[:,14]
# FAIRTEN3 = cont[:,15]
# ANCHTEN1 = cont[:,16]
# ANCHTEN2 = cont[:,17]
# ANCHTEN3 = cont[:,18]

# EV_PtfmSurge = cont_ev[:,1]
# EV_PtfmHeave = cont_ev[:,2]
# EV_PtfmPitch = cont_ev[:,3]
# EV_Wave1Elev = cont_ev[:,4]
# EV_Wave1Elev2 = cont_ev[:,5]
# EV_Wave1Elev1 = cont_ev[:,6]
# EV_HydroFxi = cont_ev[:,7]
# EV_HydroFyi = cont_ev[:,8] 
# EV_HydroFzi = cont_ev[:,9] 
# EV_HydroMxi = cont_ev[:,10] 
# EV_HydroMyi = cont_ev[:,11]  
# EV_HydroMzi = cont_ev[:,12]
# EV_FAIRTEN1 = cont_ev[:,13] 
# EV_FAIRTEN2 = cont_ev[:,14]
# EV_FAIRTEN3 = cont_ev[:,15]
# EV_ANCHTEN1 = cont_ev[:,16]
# EV_ANCHTEN2 = cont_ev[:,17]
# EV_ANCHTEN3 = cont_ev[:,18]

# plt = plot(tₑᵥ, Wave1Elev[lb:ub])
# plot!(tₑᵥ, EV_Wave1Elev)
# display(plt)

# plt = plot(tₑᵥ, PtfmSurge[lb:ub], lab = "PtfmSurge")
# plot!(tₑᵥ, EV_PtfmSurge, lab = "EV_PtfmSurge")
# display(plt)

# plt = plot(tₑᵥ, PtfmHeave[lb:ub], lab = "PtfmHeave")
# plot!(tₑᵥ, EV_PtfmHeave, lab = "EV_PtfmHeave")
# display(plt)

# plt = plot(tₑᵥ, PtfmPitch[lb:ub], lab = "PtfmPitch")
# plot!(tₑᵥ, EV_PtfmPitch, lab = "EV_PtfmPitch")
# display(plt)

# plt = plot(tₑᵥ, FAIRTEN1[lb:ub], lab = "FAIRTEN1")
# plot!(tₑᵥ, FAIRTEN2[lb:ub], lab = "FAIRTEN2")
# plot!(tₑᵥ, FAIRTEN3[lb:ub], lab = "FAIRTEN3")
# display(plt)

# plt = plot(tₑᵥ, EV_FAIRTEN1, lab = "EV_FAIRTEN1")
# plot!(tₑᵥ, EV_FAIRTEN2, lab = "EV_FAIRTEN2")
# plot!(tₑᵥ, EV_FAIRTEN3, lab = "EV_FAIRTEN3")
# display(plt)

# plt = plot(tₑᵥ, ANCHTEN1[lb:ub], lab = "ANCHTEN1")
# # plot!(tₑᵥ, ANCHTEN2[lb:ub], lab = "ANCHTEN2")
# plot!(tₑᵥ, ANCHTEN3[lb:ub], lab = "ANCHTEN3")
# display(plt)

# plt = plot(tₑᵥ, EV_ANCHTEN1, lab = "EV_ANCHTEN1")
# plot!(tₑᵥ, EV_ANCHTEN2, lab = "EV_ANCHTEN2")
# plot!(tₑᵥ, EV_ANCHTEN3, lab = "EV_ANCHTEN3")
# display(plt)

# FRev, Mev, _, _, _, Nfft, H = one_side_asp(Wave1Elev[lb:ub], tₑᵥ)
# itp_sp = interpolate(FRev, Mev, BSplineOrder(4))
# FR = range(0,fcut, Nfft)
# MAG = itp_sp.(FR)
# EV = ifft(H)

# # Save OpenFAST event
# Aₑᵥ = Wave1Elev[lb:ub]
# tₑᵥ = round.((tₑᵥ.-t[lb])*1e4)/1e4
# open("library/"*postOFpath*"/OFEV_1", "w")
# writedlm("library/"*postOFpath*"/OFEV_1", [tₑᵥ Aₑᵥ], '\t')

# plot(ηOW3D[:,1], ηOW3D[:,2])
#############################################################################################
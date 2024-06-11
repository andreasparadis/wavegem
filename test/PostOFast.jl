# Load necessary modules
using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles
using FFTW, DSP, Statistics, BSplineKit

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
pdir::String = "SE/"    # Parent directory
cdir::String = "LAF"            # Case directory
rdir::String = "/2"             # Run directory
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

# Open and read files
## Full simulation
cont = parse_fxw_pf(OFoutpath, 0, 5)
heads = readdlm("library/"*postOFpath*"/dataHdr")
## Event
cont_ev = parse_fxw_pf(postOFpath*"/outD_Ev1", 0, 5)
## Gaussian interpolation
cont_G = parse_fxw_pf(postOFpath*"/outD_Gint", 0, 5)

# ## Sum of EWGs
# cont_EWGsum = parse_fxw_pf(postOFpath*"/outData_EWGsum", 0, 5)
# ## EWG 1
# cont_EWG1 = parse_fxw_pf(postOFpath*"/outData_EWG1", 0, 5)
# ## EWG 1+2
# cont_EWG12 = parse_fxw_pf(postOFpath*"/outData_EWG1+2", 0, 5)
# ## NewWave
# cont_NW = parse_fxw_pf(postOFpath*"/outData_NW2", 0, 5)

ηOW3D = parse_fxw_pf(postOW3Dpath*"/eta", 0, 0)
ηlin = parse_fxw_pf(Dpath*"/eta_lin", 0, 0)
ηOG = parse_fxw_pf(Dpath*"/eta0", 0, 0)

#############################################################################################
Nₜ = size(cont)[1]
NoRows = size(cont)[2]

t = zeros(Float64,Nₜ)
dt = 0.0125
[t[i] = (i-1)*dt for i ∈ 1:Nₜ]

fcut = 0.625

Wave1Elev = cont[:,4]
Wave1Elev1 = cont[:,5]
# peak = findmax(PtfmSurge)
peak = findmax(Wave1Elev)
trange = Int64(round(50/dt))
lb = peak[2]-trange; ub = lb+2*trange # For event plots

# Spectra
FRof, Mof, _ = one_side_asp(Wave1Elev, t)
FRlin, Mlin, _ = one_side_asp(ηlin[:,2], ηlin[:,1])
FRog, Mog, _ = one_side_asp(ηlin[:,2], ηlin[:,1])
FR, MAG, _ = one_side_asp(Wave1Elev .- Wave1Elev1, t)

fₚ = FRof[findmax(Mof)[2]]

plt = plot(FRof, Mof, yscale=:log10, line =:dot, lab = "OpenFAST")
plot!(FRlin, Mlin, yscale=:log10, opacity=0.6, lab = "Linear")
plot!(FR, MAG, yscale=:log10, opacity=0.3, lab = "2nd order")
plot!(xlim=(0,fcut), ylim=(1e-6,1),  minorgrid=true)
display(plt)

plt = plot(FRog, Mog, lab = "Original")
plot!(FRof, Mof, line =:dash, lw = 0.5, lab = "OpenFAST")
plot!(xlim=(0, 0.625))
display(plt)

# Time histories
for i ∈ 1:12
    plti = plot(xlab = "t [s]", title = "$(heads[i])") 
    plot!(t[lb:ub], cont[lb:ub,i], lw=2, lab = "Full sim")
    plot!(t[lb:ub], cont_ev[:,i], line=:dash, lw=2, lab = "Event")
    plot!(t[lb:ub], cont_G[:,i], line=:dashdot, lw=2, lab = "Gauss WGs")
    # plot!(t[lb:ub], cont_EWGsum[:,i], lab = "sum EWGs")
    # plot!(t[lb:ub], cont_EWG1[:,i], lab = "EWG 1")
    # plot!(t[lb:ub], cont_EWG12[:,i], lab = "EWGs 1+2")
    # plot!(t[lb:ub], cont_NW[:,i], line=:dot, lw=2, lab = "NewWave")
    display(plti)
    savefig("library/"*postOFpath*"/$(heads[i]).svg")
    savefig("library/"*postOFpath*"/$(heads[i]).png")
end

# Plot RAOs
FR = readdlm("library/SE/LAF/2/0/postOFAST/FREQ", skipstart=5, '\n');
PSD = parse_fxw_pf("SE/LAF/2/0/postOFAST/PSDs", 0, 5);

pPSDs = plot(FR, PSD[:,1], xlab = L"f~[Hz]", ylab = L"PSD~[Unit^2/Hz]", lab = L"Surge", lw=3, yscale=:log10)
plot!(FR, PSD[:,2], lab = L"Heave", lw=3, yscale=:log10)
plot!(FR, PSD[:,3],  lab = L"Pitch", lw=3, yscale=:log10)
plot!([fₚ, fₚ], [1e-7, findmax(PSD[:,3])[1]], line=:dash, lw=2, color=:red, lab = L"f_p", yscale=:log10)
plot!(xlim=(0, 0.2), ylim=(1e-7, 1e3), title = "Power Spectral Densities of main DOFs", minorgrid=true)
display(pPSDs)
savefig("library/"*postOFpath*"/PSDs.svg")
savefig("library/"*postOFpath*"/PSDs.png")

T = 1 ./FR[2:end]
pPSDs = plot(T, PSD[2:end,1], xlab = L"T~[s]", ylab = L"PSD~[Unit^2/Hz]", lab = L"Surge", lw=3, yscale=:log10)
plot!(T, PSD[2:end,2], lab = L"Heave", lw=3, yscale=:log10)
plot!(T, PSD[2:end,3],  lab = L"Pitch", lw=3, yscale=:log10)
plot!([1/fₚ, 1/fₚ], [1e-7, findmax(PSD[:,3])[1]], line=:dash, lw=2, color=:red, lab = L"f_p", yscale=:log10)
plot!(xlim=(0, 80), ylim=(1e-7, 1e3), title = "Power Spectral Densities of main DOFs", minorgrid=true)
display(pPSDs)
savefig("library/"*postOFpath*"/PSDs_T.svg")
savefig("library/"*postOFpath*"/PSDs_T.png")

# PtfmSurge = cont[:,1]
# PtfmHeave = cont[:,2]
# PtfmPitch = cont[:,3]
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

# plt = plot(t[lb:ub], Wave1Elev[lb:ub])
# plot!(t[lb:ub], EV_Wave1Elev)
# display(plt)

# plt = plot(t[lb:ub], PtfmSurge[lb:ub], lab = "PtfmSurge")
# plot!(t[lb:ub], EV_PtfmSurge, lab = "EV_PtfmSurge")
# display(plt)

# plt = plot(t[lb:ub], PtfmHeave[lb:ub], lab = "PtfmHeave")
# plot!(t[lb:ub], EV_PtfmHeave, lab = "EV_PtfmHeave")
# display(plt)

# plt = plot(t[lb:ub], PtfmPitch[lb:ub], lab = "PtfmPitch")
# plot!(t[lb:ub], EV_PtfmPitch, lab = "EV_PtfmPitch")
# display(plt)

# plt = plot(t[lb:ub], FAIRTEN1[lb:ub], lab = "FAIRTEN1")
# plot!(t[lb:ub], FAIRTEN2[lb:ub], lab = "FAIRTEN2")
# plot!(t[lb:ub], FAIRTEN3[lb:ub], lab = "FAIRTEN3")
# display(plt)

# plt = plot(t[lb:ub], EV_FAIRTEN1, lab = "EV_FAIRTEN1")
# plot!(t[lb:ub], EV_FAIRTEN2, lab = "EV_FAIRTEN2")
# plot!(t[lb:ub], EV_FAIRTEN3, lab = "EV_FAIRTEN3")
# display(plt)

# plt = plot(t[lb:ub], ANCHTEN1[lb:ub], lab = "ANCHTEN1")
# # plot!(t[lb:ub], ANCHTEN2[lb:ub], lab = "ANCHTEN2")
# plot!(t[lb:ub], ANCHTEN3[lb:ub], lab = "ANCHTEN3")
# display(plt)

# plt = plot(t[lb:ub], EV_ANCHTEN1, lab = "EV_ANCHTEN1")
# plot!(t[lb:ub], EV_ANCHTEN2, lab = "EV_ANCHTEN2")
# plot!(t[lb:ub], EV_ANCHTEN3, lab = "EV_ANCHTEN3")
# display(plt)

# FRev, Mev, _, _, _, Nfft, H = one_side_asp(Wave1Elev[lb:ub], t[lb:ub])
# itp_sp = interpolate(FRev, Mev, BSplineOrder(4))
# FR = range(0,fcut, Nfft)
# MAG = itp_sp.(FR)
# EV = ifft(H)

# # Save OpenFAST event
# Aₑᵥ = Wave1Elev[lb:ub]
# tₑᵥ = round.((t[lb:ub].-t[lb])*1e4)/1e4
# open("library/"*postOFpath*"/OFEV_1", "w")
# writedlm("library/"*postOFpath*"/OFEV_1", [tₑᵥ Aₑᵥ], '\t')

# plot(ηOW3D[:,1], ηOW3D[:,2])
#############################################################################################
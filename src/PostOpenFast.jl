module PostOpenFast
using Plots, LaTeXStrings, DelimitedFiles, Dates, LinearAlgebra
using FFTW, DSP, Statistics, BSplineKit
import ColorSchemes.darkrainbow

gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))
cb = darkrainbow

include("WAVEGEM.jl")
import .WAVEGEM

# Include necessary scripts for functions
include("signal_processing.jl")
include("text_process.jl")
include("peak_detect.jl")
include("jonswap.jl")
include("directories.jl")

#############################################################################################
# Import Global Variables & Module specific inputs
ρ, g, _, d, _, _, γ, fcut = WAVEGEM.GlobInp0
pdir, run_id, phi_id, prb_id = WAVEGEM.GlobInp1
case_id, _, rundir = WAVEGEM.GlobPaths[2:4]
Decpath, _, DecFigs = WAVEGEM.GlobPaths[9:11]
runstr, case_id, casedir, rundir, phipath, OW3Dcdir, OW3Drdir, OW3Dphipath, 
            Decpath, DecEvs, DecFigs, postOFpath, postOW3Dpath, OFASTpath = WAVEGEM.GlobPaths
Wave = WAVEGEM.Wave
full, fev = WAVEGEM.POFflags
CEid = WAVEGEM.CEid
T₀ₛ¹, T₀ₛ², T₀ₕ, T₀ₚ = WAVEGEM.FOWT[:]

#############################################################################################
# Open and read files
## Full simulation
if full
    OFASTout = joinpath(postOFpath,"outD_")   # OpenFAST output file
    cont = parse_fxw_pf(OFASTout, 0, 5)
    heads = readdlm(joinpath(postOFpath,"dataHdr"))


end
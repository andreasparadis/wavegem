module Events
using Plots, LaTeXStrings, DelimitedFiles, Dates, LinearAlgebra
using FFTW, DSP, Statistics, BSplineKit

gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))

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
    
end
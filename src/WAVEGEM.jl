module WAVEGEM
# ------------------- "Frontend" of the WAVEGEM post-processing suite" -------------------- #
# --------------------- Created by Andreas Paradeisiotis, 24/11/2023 ---------------------- #

# ... To execute a module, run in the REPL, the corresponding include command below:
# >> include("src/SeaState.jl")

# export frec, flong, fplot

# Include Files Containing Necessary Functions
include("func/directories.jl")
include("func/calc_globals.jl")

#############################################################################################
# Declaration of Global Variables
## Significant wave height [m], Peak period [s], peakedness [-], Cut-off frequency [Hz]
Hₛ, Tₚ, γ, fcut::Float64 = 8.0, 10.0, 3.3, 0.625  # JONSWAP parameters  
pdir::String = "SE" # Parent directory in library folder ("SE", "WM", "rand_seas" etc.)
run_id::Int8 = 1    # Run identifier in each case (1,2,3..)
phi_id::Int8 = 0    # Phase shift for decomposition: 0, 1, 2 , 3 == rand(), -π/2, -π, -3π/2
prb_id::Int8 = 4    # ID No of probe from OW3D simulations

const ρ, g::Float64 = 1025.0, 9.81  # Sea water density [kg/m³], Accel. of gravity [m/s²]
Ldom, d::Float64 = 2000.0, 100.0    # [m] Domain length, Water depth

# Flags
flong = Bool(1) # Long simulation? Set 0: false or 1: true

# Calculated global variables
Tᵢ, λ⁻, υ⁻, tₑ = time_lims(fcut, flong)

#############################################################################################
# Declaration of Module specific inputs         

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
## ------------------------------------- SeaState -----------------------------------------##
# Module flags: 0 or 1 == false or true (Bool(0) or Bool(1))
SSflags = (Bool(1), Bool(1))    # 1: Record run?, 2: Plot results?
fₛ = 16*fcut        # Sampling frequency
Cr::Float64 = 0.625 # Courant number
Aid::Int8 = 0      # Amplitude coeffs - 0:ones(), 1:Normal, 2:Rayleigh, Other:abs(real(amp))

SeaStateInp = (fₛ, Cr, Aid)

## ----------------------------------- Decomposition --------------------------------------##
# Module flags: 0 or 1 == false or true (Bool(0) or Bool(1))
Dflags = (Bool(1), Bool(1))    # 1: Record run?, 2: Plot results?
# Signal truncation
tstart = round(1000/υ⁻)  # Depending on position of FOWT (x=1000m) and slowest wave [s]
tend = tₑ                # Adjust accordingly [s]
t_range = round(sqrt(2π/g * Ldom/2)) # Event length [s]
NP = 5  # Number of peaks to process (== No of events)

DecSigs = ("eta0", "eta1", "eta2", "eta3") # Signal files for decomposition (~/Decomposition/)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# For import to other modules
GlobInp0 = (ρ, g, Ldom, d, Hₛ, Tₚ, γ, fcut)
GlobInp1 = (pdir, run_id, phi_id, prb_id)
GlobPaths = paths(pdir, run_id, phi_id, prb_id, flong)  # Paths

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

# push!(LOAD_PATH, "/src/")
end
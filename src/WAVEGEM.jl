module WAVEGEM
# ------------------- "Frontend" of the WAVEGEM post-processing suite" -------------------- #
# --------------------- Created by Andreas Paradeisiotis, 24/11/2023 ---------------------- #

# ... To execute a module, run in the REPL, the corresponding include command below:
# >> include("src/SeaState.jl")

# export frec, flong, fplot

# Include Files Containing Necessary Functions
include("func/directories.jl")

#############################################################################################
## Declaration of Global Variables
# Significant wave height [m], Peak period [s], peakedness [-], Cut-off frequency [Hz]
Hₛ, Tₚ, γ, fcut::Float64 = 3.0, 8.0, 3.3, 0.625  # JONSWAP parameters  
pdir::String = "SE" # Parent directory in library folder ("SE", "WM", "rand_seas" etc.)
run_id::Int8 = 1    # Run identifier in each case (1,2,3..)
phi_id::Int8 = 0    # Phase shift for decomposition: 0, 1, 2 , 3 == rand(), -π/2, -π, -3π/2
prb_id::Int8 = 5    # ID No of probe from OW3D simulations

const ρ, g::Float64 = 1025.0, 9.81  # Sea water density [kg/m³], Accel. of gravity [m/s²]
Ldom, d::Float64 = 2000.0, 100.0    # [m] Domain length, Water depth

# Flags
flong = Bool(0) # Long simulation? Set 0: false or 1: true

#############################################################################################
# Calculated global variables
Tᵢ = 1/fcut             # Cut-off period of JONSWAP spectrum [s]
λ⁻ = g/(2*π) * Tᵢ^2     # Shortest wavelength (Deep water) [m]
υ⁻ = λ⁻/Tᵢ              # Slowest wave [m/s]

if flong 
    tₑ = 180*60 # 3-hour duration of simulation [s]
else
    tₑ = round(2*Ldom/υ⁻ /10)*10 # Duration of simulation [s]
end

#############################################################################################
# Declaration of Module specific inputs

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
## ------------------------------------- SeaState -----------------------------------------##
fₛ = 16*fcut        # Sampling frequency
Cr::Float64 = 0.625 # Courant number
A_id::Int8 = 0      # Amplitude coeffs - 0:ones(), 1:Normal, 2:Rayleigh, Other:abs(real(amp))

SSflags = (Bool(0), Bool(1))    # 1: Record run?, 2: Plot results?, Set 0: false or 1: true
SeaStateInp = (fₛ, Cr, A_id)

## ----------------------------------- Decomposition --------------------------------------##
# Signal truncation
tstart = round(1000/υ⁻)  # Depending on position of FOWT (x=1000m) and slowest wave [s]
tend = tₑ                # Adjust accordingly [s]

DecSigs = ("eta0", "eta1", "eta2", "eta3") # Signal files for decomposition (~/Decomposition/)
Dflags = (Bool(1), Bool(1))    # 1: Record run?, 2: Plot results?, Set 0: false or 1: true

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
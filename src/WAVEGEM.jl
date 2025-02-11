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
Hₛ, Tₚ, γ::Float64 = 5.0, 8.0, 3.3  # JONSWAP parameters
fₛ, Tₑ::Float64 = 2^3, 2^5           # Sampling frequency and return period 
pdir::String = "SE" # Parent directory in library folder ("SE", "WM", "rand_seas" etc.)
run_id::Int8 = 15    # Run identifier in each case (1,2,3..)
phi_id::Int8 = 0    # Phase shift for decomposition: 0, 1, 2 , 3 == rand(), -π/2, -π, -3π/2
prb_id::Int8 = 4    # ID No of probe from OW3D simulations

const ρ, g::Float64 = 1025.0, 9.81  # Sea water density [kg/m³], Accel. of gravity [m/s²]
Ldom, d::Float64 = 2000.0, 100.0    # [m] Domain length, Water depth

# Flags
flong = true # Long simulation? Set true or false
# Global variables
Tᵢ, λ⁻, υ⁻, tₑ, fcut = time_lims(flong, Tₑ)
GlobInp0 = (ρ, g, Ldom, d, Hₛ, Tₚ, γ, fcut)
GlobInp1 = (pdir, run_id, phi_id, prb_id)
GlobPaths = paths(pdir, run_id, phi_id, prb_id, flong)  # Paths

#############################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODULE SPECIFIC INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#         
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
## ------------------------------------- SeaState -----------------------------------------##
# Module flags: true or false 
SSflags = (true, true)    # 1: Record run?, 2: Plot results?
Nₛ = Int(tₑ/2)
Cr::Float64 = 0.625       # Courant number
Aid::Int8 = 0      # Amplitude coeffs - 0:ones(), 1:Normal, 2:Rayleigh, Other:abs(real(amp))

SeaStateInp = (fₛ, Nₛ, Cr, Aid)
## ----------------------------------------------------------------------------------------##

## ----------------------------------- Decomposition --------------------------------------##
# Module flags: true or false 
Dflags = (true, true)    # 1: Record run?, 2: Plot results?
# Signal truncation
tstart = round(1000/υ⁻)  # Depending on position of FOWT (x=1000m) and slowest wave [s]
tend = tₑ-0.25                # Adjust accordingly [s]
t_range = round(sqrt(2π/g * Ldom/2)) # Event length [s]
NP = 5  # Number of peaks to process (== No of events)

DecSigs = ("eta0", "eta1", "eta2", "eta3") # Signal files for decomposition (~/Decomposition/)

# -------------------------------------- Events ------------------------------------------- #
# Module flags: true or false 
Evflags = (true, true, true)    # 1: Record?, 2: Plot?, 3: Suppress Decomposition rec & plot?
# Suppress record & plot for Decomposition module if it has already been applied
fid_1st = joinpath(GlobPaths[1],"eta_lin")
if Evflags[3] && !isfile(fid_1st)
    Dflags = (false,false)
end

## ----------------------------------------------------------------------------------------##

# ---------------------------------- PostOpenFast ----------------------------------------- #
# Module flags: true or false 
POFflags = (true,false)     # 1: Full simulation?, 2: Partial (event) simulation?
CEid = 0  # Critical event type = 0: Max Fair. tension, 1: Max Pitch, 2: Max Wave, 3: Max CoM 

## ----------------------------------------------------------------------------------------##


## -------------------------------------- ERGEEs ------------------------------------------##
evdir = "MaxFair_10"    # Event directory
ERGEEflags = (true, false)  # 1: Record run?, 2: Truncate signal?

# Trunc.: false ≡ 1st to last upcross., true ≡ Specific upcross. prior and after highest peak
pNuc = 4;   aNuc = 2    # Relevant only for true
# Envelope low-pass filter
fcˡᵖ = 0.5*fcut # Cut-off frequency
# Envelope peak detection
std_fac = 0 # MinPeakVal - std coeff (1.25: >η̅, 2: >Hₛ/2, 0: all)
MinPeakDist = 2*Tᵢ
# Runge-Kutta solution parameters
RK_ϵ = 1e-5        # Accuracy
RK_Nτ = Int(1e4)   # Max iterations
RK_dτ = 1e-1       # Step of fictitious time
RK_Tlb = 0*Tᵢ    # Constraint - Lower bound 
RK_Tub = Tₑ # Constraint - Higher bound

ERGEEsInp = (fcˡᵖ, std_fac, MinPeakDist, RK_ϵ, RK_Nτ, RK_dτ, RK_Tlb, RK_Tub)
## ----------------------------------------------------------------------------------------##

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
## FOWT eigenfrequencies
whichFOWT = 2   # 1: OC4 semi-submersible, 2: Mikes
FOWT = [0;0;0;0]    # T₀ₛ¹, T₀ₛ², T₀ₕ, T₀ₚ 
if whichFOWT == 1
    FOWT[:] = [250.0000; 35.714107; 17.543772; 27.026892] 
elseif whichFOWT == 2
    FOWT[:] = [93.079545; 26.944079; 17.064583; 26.944079]
end

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
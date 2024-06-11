using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles, FFTW
using Statistics, Random, Distributions
using Dates, SparseArrays

include("wave_theory.jl")
include("signal_processing.jl")
include("jonswap.jl")
include("peak_detect.jl")
include("directories.jl")

#############################################################################################
# Input variables
Hₛ, Tₚ, γ::Float64 = 1.0, 8.0, 3.3 # JONSWAP parameters
Ldom, d::Float64 = 2000.0, 100.0 # [m] Domain length, Water depth
ρ, g::Float64 = 1025.0, 9.81

mutable struct Wave 
    f::Float64
    T::Float64
    ω::Float64
    λ::Float64
    κ::Float64
    υᶜ::Float64
end

## Simulation details
frec = Bool(0)  # Record run?      Set 0: false or 1: true
flong = Bool(0) # Long simulation? Set 0: false or 1: true
fplot = Bool(0) # Plot results?    Set 0: false or 1: true

fcut = 0.625  # [Hz] Cut-off frequency
fₛ = 16*fcut  # Sampling frequency
Cr = 0.625    # Courant number
A_id = 0      # Amplitude coeffs - 0: ones(), 1: Normal, 2: Rayleigh, Other: abs(real(amp))
#############################################################################################

if frec
    pdir::String = "SE" # Directory name in library folder (WM, SE, rand_seas etc.)
    run_id::Int8 = 2    # Run identifier in each case (1-10)
    phi_id::Int8 = 0    # Phase shift for decomposition: 0, 1, 2 , 3 == rand(), -π/2, -π, -3π/2

    runstr, rundir, phipath, OW3Drdir, OW3Dinp, _ = paths(pdir, run_id, phi_id, 5, flong)

    ## Info for creating directories and writing variables in text files
    head_o3d = ["# JS_$(Hₛ)m_$(Tₚ)s_ϕ$phi_id" " "] # Output files header

    # Print useful information
    println("Paths:")
    println("Run folder: $runstr")
    println("Path to OW3D input folder: $OW3Dinp")
    println("Case path in library folder: $phipath")

    # Make new directories
    if isdir(phipath)
        error("Andreas: This folder already exists, since this run has already been recorded!")
    elseif isdir(rundir)
        mkdir(phipath)
        mkdir(OW3Dinp)
    else
        mkdir(rundir)
        mkdir(phipath)
        mkdir(OW3Drdir)
        mkdir(OW3Dinp)
    end
else
    phi_id::Int8 = 0    # Random phase
    println("Run is not recorded.")
end

#############################################################################################
## Definition of temporal and spatial limits
ωcut = 2π*fcut # [rad/s]
Tᵢ = 1/fcut     # Cut-off frequency for JONSWAP spectrum

# Wave characteristic values at peak period Tp
peak_wave = wave_qnts(Tₚ,d)
fₚ = peak_wave.f;   ωₚ = peak_wave.ω  
λₚ = peak_wave.λ;   κₚ = peak_wave.κ
υₚ = peak_wave.υᶜ

# Determine maximum period (deep water assumption does not apply)
Tₑ₀ = sqrt(2π/g * Ldom/2) # Initial guess for maximum period

long_wave = wave_qnts(Tₑ₀,d)
f⁻ = long_wave.f;   ω⁻ = long_wave.ω
λ⁺ = long_wave.λ;   κ⁻ = long_wave.κ
υ⁺ = long_wave.υᶜ

Tₑ = round(2π/ω⁻) # Maximum period

# Duration of simulation - Slowest wave returns to origin
short_wave = wave_qnts(Tᵢ,d)
λ⁻ = short_wave.λ # Minimum wavelength
ω⁺ = short_wave.ω    
κ⁺ = short_wave.κ    
υ⁻ = short_wave.υᶜ

if flong 
    tₑ = 180*60 # 3-hour simulation
else
    tₑ = round(2*Ldom/υ⁻ /10)*10
end

# Sampling frequency fs>=2*fcut=fnyq
fnyq = 2*fcut
ωₛ = 2π*fₛ
dt = 1/fₛ
Nₜ = Int64(round(tₑ/dt))+1
# N = 2*ncut-1

# Discretization of JONSWAP frequency range
max_kd = κ⁺*d
df = 1/((Nₜ-1)*dt); dω = 2π*df  
Ncut = Int64(round((ω⁺-ω⁻)/dω))
# Nₛ = Int64(round(Nₜ/2))+1
Nₛ = Ncut

fⱼ,Sⱼ = spectrum(Hₛ,Tₚ,γ,Tᵢ,Tₑ,Nₛ) # Generate JONSWAP spectrum for pair (Hₛ,Tₚ) ∈ [Tᵢ,Tₑ]
η̂, t, dt, ηₜ, u, u̇, ẇ, Φ, ηₜᵈ, ϕ, x, dx, Nₓ, ηₓ, z, dz, Nz, U = rand_sea(fⱼ, Sⱼ, Tₑ, dt, Nₜ, phi_id, A_id, Ldom, d, Cr)

#############################################################################################
# Simulation info
if ωcut*dt >= 2*sqrt(2)
    print("OCW3D WARNING: The Runge-Kutta stability condition is not satisfied: ωcut•dt=",ωcut*dt ," ≥ 2√2 \n")
else
    print("OCW3D: The Runge-Kutta stability condition is satisfied: ωcut•dt=",ωcut*dt ," ≤ 2√2 \n")
end

content = ["JONSWAP: Hₛ=$Hₛ [m] , Tₚ=$Tₚ [sec] , γ=$γ , d=$d [m]";
        "T ∈ $([Tᵢ Tₑ]) [sec]";
        "f ∈ $([round(f⁻*1e4)/1e4 fcut]) [Hz]";
        "ω ∈ $([round(ω⁻*1e4)/1e4 round(ωcut*1e4)/1e4]) [rad/s]";
        "λ ∈ $([round(λ⁻*1e4)/1e4 round(λ⁺*1e4)/1e4]) [m]";
        "κ ∈ $([round(κ⁻*1e4)/1e4 round(κ⁺*1e4)/1e4]) [1/m]";
        "κd ∈ $([κ⁻*d max_kd]) [-]";
        "dt = $dt [s]";
        "dx = $dx [m]";
        "dz = $dz [m]";
        "Duration: $((Nₜ-1)*dt) [sec]";
        "Nₜ = $Nₜ";
        "Nₓ = $Nₓ";
        "Nₛ = $Nₛ";
        "Nz = $Nz";
        "Number of spatial points for smallest wavelength: $(λ⁻/dx)";
        "Hₛ of generated sea-state: $(round(4*std(ηₜ) *1e4)/1e4) [m]"]

lcont = length(content)
for i ∈ 1:lcont
    println(content[i])
end

if frec
    fid::String = "INFO" # File name
    open(phipath*fid, "w")
    writedlm(phipath*fid, content, '\t')

    info1 = [string(dt)*" "*string(Nₜ)*" 1" " "]
    info2 = [string(1) " "]

    # Surface elevation (temporal)
    fid::String = "eta_t" # File name
    open(phipath*fid, "w")
    head = ["t [s]" "η [m]"]
    writedlm(phipath*fid, [head; t ηₜ], '\t')

    # Surface elevation for OW3D wave file input (no time column)
    fid::String = "eta" # File name
    open(OW3Dinp*fid, "w")
    writedlm(OW3Dinp*fid, [head_o3d[1]; dt; ηₜ], '\t')

    # Phase
    fid::String = "phi" # File name
    open(phipath*fid, "w")
    head = ["f [Hz]" "ϕ [-]"]
    writedlm(phipath*fid, [head; fⱼ ϕ], '\t')

    # Horizontal velocity
    fid::String = "u" # File name
    head = ["t [s]" "u [m/s]"]
    open(phipath*fid, "w")
    writedlm(phipath*fid, [head; t u], '\t')
    # For OW3D wavemaker signal (flux)
    open(OW3Dinp*fid, "w")
    writedlm(OW3Dinp*fid, [head_o3d; info1; info2; t u], '\t')

    # Potential
    fid::String = "Potential" # File name
    open(phipath*fid, "w")
    writedlm(phipath*fid, [head_o3d; info1; info2; t Φ], '\t')

    # # Uniform horizontal velocity
    # fid::String = "U" # File name
    # open(phipath*fid, "w")
    # writedlm(phipath*fid, [head_o3d; info1; info2; t U], '\t')

    # # Surface elevation (spatial)
    # fid::String = "eta_x" # File name
    # open(phipath*fid, "w")
    # head = ["x [m]" "η [m]"]
    # writedlm(phipath*fid, [head; x ηₓ], '\t')

    # # Flux
    # fid::String = "Flux" # File name
    # open(phipath*fid, "w")
    # writedlm(phipath*fid, [head_o3d; info1; info2; t flux], '\t')

    # # Normalized surface elevation (temporal)
    # fid::String = "eta_t_norm" # File name
    # open(phipath*fid, "w")
    # writedlm(phipath*fid, [head_o3d; info1; info2; t ηₜ/(findmax(ηₜ)[1])], '\t')

    # # Horizontal velocity at z = -1.5 m
    # fid::String = "ud" # File name
    # open(phipath*fid, "w")
    # writedlm(phipath*fid, [head_o3d; info1; info2; t ud], '\t')

    # # Wavemaker displacement
    # fid::String = "WM_z" # File name
    # open(phipath*fid, "w")
    # writedlm(phipath*fid, [head_o3d; info1; info2; t ζ], '\t')
end

#############################################################################################
# PLOTS
if fplot
    plt = plot(t, ηₜ, xlab = L"t~[s]", ylab = L"\eta~[m]", lab = L"η(t)", lw=1)
    display(plt)
    if frec
        savefig(phipath*"eta_t.svg")
    end

    # plt = plot(t, ηₜᵈ, xlab = L"t~[s]", ylab = L"\eta (z=-d,t)~[m]", lab = L"η(t)", lw=1)
    # display(plt)

    plt = plot(t, ηₜᵈ/(findmax(ηₜᵈ)[1]), lab = L"\overline{\eta}(z=-d)", lw=2)
    plot!(t, ηₜ/(findmax(ηₜ)[1]), xlab = L"t~[s]", ylab = L"\overline{\eta}~[-]", lab = L"\overline{\eta}(z=0)", lw=1)
    # plot!(xlim=(500,1000))
    display(plt)

    plt = plot(t, u, xlab = L"t~[s]", ylab = L"u~[m/s]", lab = L"u(t)", lw=1) 
    display(plt)

    plt = plot(t, ẇ, xlab = L"t~[s]", ylab = L"ẇ~[m/s]", lab = L"ẇ(t)", lw=1) 
    plot!([t[1]; t[end]], [0.5*g; 0.5*g], lab = L"Breaking~limit", line=:dash)
    display(plt)

    # plt = plot(t, U, xlab = L"t~[s]", ylab = L"U~[m/s]", lab = L"U(t)", lw=1) 
    # display(plt)

    # plt = plot(t, flux, xlab = L"t~[s]", ylab = L"Flux~[m^2/s]", lab = L"u(z=-d)~d", lw=1) 
    # display(plt)

    # plt = plot(t, vflux, xlab = L"t~[s]", ylab = L"Volume Flux~[m^3/s]", lab = L"u(z=-d)~d", lw=1) 
    # display(plt)

    # plt = plot(t, P, xlab = L"t~[s]", ylab = L"\overline{P}~[W/m]", lab = L"η(t)", lw=2) 
    # display(plt)

    freq, mag, _, dfₜ = one_side_asp(ηₜ,t) # Surface elevation spectrum
    plt = plot(freq, mag, xlab = L"f~[Hz]", ylab = L"\tilde{\eta}~[m]", lab = "Non-linear spectrum", lw=1)
    plot!(fⱼ, η̂, xlab = L"f~[Hz]", ylab = L"\tilde{\eta}~[m]", lab = "JONSWAP component amplitudes", lw=1)
    plot!(xlim=(0,fⱼ[end]))
    display(plt)

    plt = plot(freq, mag.^2 /(2*dfₜ), xlab = L"f~[Hz]", ylab = "Magn", lab = "Non-linear variance spectrum", lw=1)
    plot!(fⱼ, Sⱼ, xlab = L"f~[Hz]", ylab = L"S(f)~[m^2/Hz]", lab = "JONSWAP variance spectrum", lw=1)
    plot!(xlim=(0,fⱼ[end]))
    display(plt)

    # freq, mag, _, _ = one_side_asp(U,t)
    # plot(freq, mag, xlab = L"f~[Hz]", ylab = L"\tilde{U}~[m]", lab = "Spectrum U", lw=1)
    # plot!(fⱼ, ω.*η̂, xlab = L"f~[Hz]", ylab = L"\tilde{\eta}~[m]", lab = "JONSWAP component amplitudes", lw=1)

    # plt = plot(x, ηₓ, xlab = L"x~[m]", ylab = L"\eta~[m]", lab = L"η(x)", lw=1, aspect_ratio=80)
    # display(plt)
    # savefig("eta_x.svg")

    # wnum, magₓ, _, dfₓ = one_side_asp(real(ηₓ),x)
    # plot(2π*wnum, magₓ, xlab = L"k~[rad/m]", ylab = L"\tilde{\eta}~[m^2]", lab = "Spatial spectrum", lw=1)
end
#############################################################################################
# Deallocate large matrices/vectors
t, ηₜ, u, u̇, ẇ, Φ = (nothing for _ = 1:6)
GC.gc()

#############################################################################################
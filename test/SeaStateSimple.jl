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
Hₛ, Tₚ, γ::Float64 = 3.0, 8.0, 3.3 # JONSWAP parameters
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
frec = false  # Record run?      Set 0: false or 1: true
fplot = true # Plot results?    Set 0: false or 1: true

fₛ = 2^3            # [Hz] Sampling frequency
Tₑ = 2^6            # Maximum period / Return period
tₑ = 2^7*Tₑ         # [s] Simulation duration
df = 1/tₑ           # Frequency resolution
Nₛ = Int(tₑ/2)      # No of frequency components
fcut = Nₛ*df + 1/Tₑ # [Hz] Cut-off frequency (1/Tₑ+1/2)

A_id = 0      # Amplitude coeffs - 0: ones(), 1: Normal, 2: Rayleigh, Other: abs(real(amp))

#############################################################################################

if frec
    phi_id::Int8 = 0    # Phase shift for decomposition: 0, 1, 2 , 3 == rand(), -π/2, -π, -3π/2

    ## Info for creating directories and writing variables in text files
    ProjPath::String = pwd()               # Project path
    libpath = joinpath(ProjPath,"library") # Path to project's library folder
    pdir::String = joinpath(libpath,"UCL","Towing Tank Tests","AP_FR") # Directory name in library folder
    
    phi_tag = ("_00.fronts", "_05.fronts", "_10.fronts", "_15.fronts")
    suffix = "HS00$(Int(round(100*Hₛ)))TP$(Int(round(Tₚ)))TR$Tₑ" # Suffix of output files
    ffront = suffix*phi_tag[phi_id+1]

    # Print useful information
    println("Results stored in:")
    println(pdir)

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
long_wave = wave_qnts(Tₑ,d)
f⁻ = long_wave.f;   ω⁻ = long_wave.ω
λ⁺ = long_wave.λ;   κ⁻ = long_wave.κ
υ⁺ = long_wave.υᶜ

# Duration of simulation - Slowest wave returns to origin
short_wave = wave_qnts(Tᵢ,d)
λ⁻ = short_wave.λ # Minimum wavelength
ω⁺ = short_wave.ω    
κ⁺ = short_wave.κ    
υ⁻ = short_wave.υᶜ

# Sampling frequency fs>=2*fcut=fnyq
fnyq = 2*fcut
ωₛ = 2π*fₛ
dt = 1/fₛ
Nₜ = Int64(round(tₑ/dt))

# Discretization of JONSWAP frequency range
max_kd = κ⁺*d
# df = 1/((Nₜ-1)*dt); dω = 2π*df  
# Nₛ = Int64(round((ω⁺-ω⁻)/dω))+1

fⱼ,Sⱼ = spectrum(Hₛ,Tₚ,γ,Tᵢ,Tₑ,Nₛ) # Generate JONSWAP spectrum for pair (Hₛ,Tₚ) ∈ [Tᵢ,Tₑ]

ω = 2π * fⱼ;    k = 1/g * ω.^2;     #υ = ω./k; υ[1]=0
sp_peak = findmax(Sⱼ);     iₚ = sp_peak[2];    fₚ = fⱼ[iₚ]

for i ∈ 1:iₚ
    if 2π/ω[i] < Tₑ
        λ,_ = dispersion(2π/ω[i],d)
        k[i] = 2π/λ
    end
end

# Initialize vectors
t = zeros(Float32,Nₜ)
[t[i] = (i-1)*dt for i ∈ 2:Nₜ]

𝚶⃗ₜ = zeros(Float32,Nₜ)

ηₜ, u, u̇, ẇ, Φ = (𝚶⃗ₜ[:] for _ = 1:5)
# Wave component amplitudes
η̂ = sqrt.(2*Sⱼ*df) 
m₀ = Hₛ^2/16
H₀ = sqrt(2*m₀*log(Nₛ))

# Phase
amp = rand(Normal(),Nₛ) + 1im*rand(Normal(),Nₛ) # Random complex vector (Normal distribution, 0 mean, 1 variance)
if phi_id == 0
    ϕ = angle.(amp) # ϕ ∈ [0,2π]
else
    global ϕ .+= π/2
    print("-$(phi_id)π/2 phase shift \n")
end

# Randomly distributed amplitude coefficient
if A_id == 0
    A = ones(Float64,Nₛ)
    println("Unit amplitude coefficients")
elseif A_id == 1
    A = rand(Normal(),Nₛ)
    println("Normally distributed amplitude coefficients")
elseif A_id == 2
    seed = Int64(round(rand(1)[1]*1000));   Random.seed!(seed);    
    distr = Rayleigh(1);    A = rand(distr,Nₛ) # Rayleigh distribution with unit scale
    println("Rayleigh distributed amplitude coefficients")
else
    A = abs.(real(amp))
    println("Amplitude coefficients from normally distributed real part of amp")
end

# Surface Wave kinematics at x=0
Ĥ = A .* η̂
Û = ω .* Ĥ
α̂ = ω.^2 .* Ĥ
Φ̂ = g * Ĥ ./ ω

for m ∈ 1:Nₛ
    θ = -ω[m]*t .+ ϕ[m]

    csθ = cos.(θ)
    snθ = sin.(θ)

    ηₜ[:] += Ĥ[m] * csθ

    Φ[:] += Φ̂[m] *snθ
    u[:] += Û[m] * csθ
    u̇[:] += α̂[m] * snθ
    ẇ[:] += - α̂[m] * snθ
end

#############################################################################################
# Simulation info
content = ["JONSWAP: Hₛ=$Hₛ [m] , Tₚ=$Tₚ [sec] , γ=$γ , d=$d [m]";
        "T ∈ $([Tᵢ Tₑ]) [sec]";
        "f ∈ $([round(f⁻*1e4)/1e4 fcut]) [Hz]";
        "ω ∈ $([round(ω⁻*1e4)/1e4 round(ωcut*1e4)/1e4]) [rad/s]";
        "λ ∈ $([round(λ⁻*1e4)/1e4 round(λ⁺*1e4)/1e4]) [m]";
        "κ ∈ $([round(κ⁻*1e4)/1e4 round(κ⁺*1e4)/1e4]) [1/m]";
        "κd ∈ $([κ⁻*d max_kd]) [-]";
        "dt = $dt [s]";
        "Duration: $((Nₜ-1)*dt) [sec]";
        "Nₜ = $Nₜ";
        "Nₛ = $Nₛ";
        "Hₛ of generated sea-state: $(round(4*std(ηₜ) *1e4)/1e4) [m]"]

lcont = length(content)
for i ∈ 1:lcont
    println(content[i])
end

#############################################################################################
# PLOTS
if fplot
    plt = plot(t, ηₜ, xlab = L"t~[s]", ylab = L"\eta~[m]", lab = L"η(t)", lw=1)
    display(plt)
    if frec
        savefig(joinpath(pdir,"eta_t.png"))
    end

    plt = plot(t, u, xlab = L"t~[s]", ylab = L"u~[m/s]", lab = L"u(t)", lw=1) 
    display(plt)

    plt = plot(t, ẇ, xlab = L"t~[s]", ylab = L"ẇ~[m/s]", lab = L"ẇ(t)", lw=1) 
    plot!([t[1]; t[end]], [0.5*g; 0.5*g], lab = L"Breaking~limit", line=:dash)
    display(plt)

    freq, mag, ang, dfₜ = one_side_asp(ηₜ,t) # Surface elevation spectrum

    plt = plot(freq, mag, xlab = L"f~[Hz]", ylab = L"\tilde{\eta}~[m]", lab = "Non-linear spectrum", lw=1)
    plot!(fⱼ, η̂, xlab = L"f~[Hz]", ylab = L"\tilde{\eta}~[m]", lab = "JONSWAP component amplitudes", lw=1)
    plot!(xlim=(0,fⱼ[end]))
    display(plt)
    if frec
        savefig(joinpath(pdir,"amp_spec.png"))
    end

    plt = plot(freq, mag.^2 /(2*dfₜ), xlab = L"f~[Hz]", ylab = "Magn", lab = "Non-linear variance spectrum", lw=1)
    plot!(fⱼ, Sⱼ, xlab = L"f~[Hz]", ylab = L"S(f)~[m^2/Hz]", lab = "JONSWAP variance spectrum", lw=1)
    plot!(xlim=(0,fⱼ[end]))
    display(plt)
    if frec
        savefig(joinpath(pdir,"var_spec.png"))
    end
end

#############################################################################################
# OUTPUT
if frec
    fid = joinpath(pdir,"INFO.txt") # File name
    open(fid, "w")
    writedlm(fid, content, '\t')

    # .fronts file
    fid = joinpath(pdir,ffront) 
    open(fid, "w")
    head0 = [suffix "" "" ""]
    head1 = ["f" "a" "angle" "phase"]
    head2 = ["Hz" "m" "rad" "rad"]
    wcont = [head0; head1; head2; fⱼ η̂ zeros(Float64,Nₛ) ϕ]
    # wcont = [head0; head1; head2;freq mag zeros(Float64,length(freq)) ang]
    writedlm(fid, wcont, '\t')

    # Surface elevation (temporal)
    fid = joinpath(pdir,"ETA_"*suffix*".txt") 
    open(fid, "w")
    head0 = [suffix ""]
    head1 = ["Time" "Height"]
    head2 = ["s" "m"]
    wcont2 = [head0; head1; head2; t ηₜ]
    writedlm(fid, wcont2, '\t')

    # # Phase
    # fid::String = "phi" # File name
    # open(pdir*fid, "w")
    # head = ["f [Hz]" "ϕ [-]"]
    # writedlm(pdir*fid, [head; fⱼ ϕ], '\t')
end
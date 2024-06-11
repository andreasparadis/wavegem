## Post-processing of the output files of OW3D simulations
# AP: This works for wave generation using the JONSWAP spectrum. 
# AP: In other cases might require some tweaks depending on the format of output files.

## Load necessary modules
using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles
using FFTW, DSP, Dates, Statistics
using CurveFit: curve_fit, RationalPoly

## Include necessary scripts for functions
include("text_process.jl")
include("signal_processing.jl")
include("jonswap.jl")
include("peak_detect.jl")
include("directories.jl")

## Read necessary text files in library folder
# Simulation Identifier
sim_id::String = "WM_run01" # Corresponding folder in library
# Definition of paths to other files
inp_fid::String = sim_id*"/OceanWave3D.inp" # Path to input file
elev_fid::String = sim_id*"/eta15" # Path to wave elevation file
spec_fid::String = sim_id*"/spectrum"       # Path to wave generation spectrum file
recon_fid::String = sim_id*"/eta0_coeffs"   # Path to signal reconstruction file

Lₓ, Nₓ, _, Tₚ, Hₛ, γ, d, max_kd = read_ow3d_inp(inp_fid)
elevation = parse_fxw_pf(elev_fid, 0, 0)    # [t η]
# open("library/"*elev_fid, "r")   
# elevation = readdlm("library/"*elev_fid, ' ', Float64, skipstart=1, '\n')
WGspec = parse_fxw_pf(spec_fid, 0, 1)       # [f T S]

# Original time vector and surface elevation
t, η, N = elevation[:,1], elevation[:,2], length(elevation[:,1])
dt = (t[end]-t[1])/(N-1)
H₁₋₃ = 4*std(η)      

# Wave generation spectrum in OW3D
f_ow3, T_ow3, S_ow3, N_ow3 = WGspec[:,1], WGspec[:,2], WGspec[:,3], length(WGspec[:,1])
spec_plt = plot(f_ow3, S_ow3, xlab = L"f~[Hz]", ylab = L"S(f)~[m^2 s]", lab = L"OW3D", lw=3)

## Generate JONSWAP spectrum for pair (Hₛ,Tₚ) ∈ [Tᵢ,Tₑ]
Tᵢ, Tₑ = 0.4, 1000.0    # Max and min period for spectrum
fᵢ = 1/Tᵢ; fₑ = 1/Tₑ    # Min and cut-off frequency
λ_max = 2π / (max_kd/d) # Max wavelength 
fⱼ, Sⱼ = spectrum(Hₛ,Tₚ,γ,Tᵢ,Tₑ,1000)
plot!(fⱼ, Sⱼ, xlab = L"f~[Hz]", ylab = L"S(f)~[m^2 s]", lab = L"JONSWAP", lw=3)
scat = height_scatter(0,25,100) # Give peak period range for the scatter diagram

## Non-linear spectrum - Entire signal
fₙₗₙ, Aₙₗₙ, _, dfₙₗₙ = one_side_asp!(η,t)
fₙₗₙ_id = floor(Integer,1/dfₙₗₙ)

## Truncated signal for processing 
Mₛ = Int(floor(540/dt));
Mₑ = Int(floor(640/dt));
# Assign truncated signal to new vars (indice 'm' from length)
tₘ = t[Mₛ:Mₑ];  ηₘ = η[Mₛ:Mₑ];  m = length(ηₘ)
# Save truncated signal to new .txt
fid_new = sim_id*"_trunc_$Mₛ-$Mₑ.txt" # Path to input file
open("library/"*fid_new, "w")
writedlm("library/"*fid_new, [tₘ ηₘ], '\t')

# Spectrum of truncated signal
freq, mag, _, df = one_side_asp!(ηₘ,tₘ)
# Spectrum plotted up to specified frequency
# Truncation at 0.5 Hz
f_index = floor(Integer,1/df)
lbl_x = L"f~[Hz]"; lbl_y = "Magnitude"; llbls = "Original signal"

nlnspec_plt = plot(fₙₗₙ[1:fₙₗₙ_id], Aₙₗₙ[1:fₙₗₙ_id], xlab=lbl_x, ylab=lbl_y, lab = llbls, lw=3)
# nlnspec_plt = plot(fₙₗₙ[1:fₙₗₙ_id], Aₙₗₙ[1:fₙₗₙ_id].^2 /(2*dfₙₗₙ), xlab=lbl_x, ylab=lbl_y, lab = llbls, lw=3)
plot!(freq[1:f_index], mag[1:f_index], xlab=lbl_x, ylab=lbl_y, lab="Truncated signal", lw=1, line=:dash)

# Hilbert transform and envelope of surface elevation
ηᴴ = hilbert(ηₘ);   ηᴴₑₙᵥ = abs.(ηᴴ)

# Peak detection based on 2nd derivative of signal
tₘ⁺, tₘ¯, ηₘ⁺, ηₘ¯, η̇ₘ, η̈ₘ = peaks(ηₘ,tₘ)
# tₘ⁺, tₘ¯, ηₘ⁺, ηₘ¯, η̇ₘ, η̈ₘ = peaks(ηᴴₑₙᵥ,tₘ)

# # Nonlinear spectrum fitting
# fit = curve_fit(RationalPoly, freq[1:f_index], mag[1:f_index], 1, 3)
# sp_fit = zeros(Float64, f_index)
# for i ∈ 1:f_index
#     sp_fit[i] = fit(freq[i])
# end

#############################################################################################
# PLOTS

# Truncated signal and envelope
lbl_x = L"t~[s]";   lbl_y = L"η~[m]";   llbls = [L"η (t)" L"Envelope"]
plt = plot(tₘ, [ηₘ ηᴴₑₙᵥ], xlab=lbl_x, ylab=lbl_y, lab=llbls, lw=[1 2])
display(plt)

# Scatter plot of extrema
plot(tₘ, ηₘ, xlab = L"t~[s]", ylab = L"\eta~[m]", lab = L"η(t)")
# plot(tₘ, ηᴴₑₙᵥ, xlab = L"t~[s]", ylab = L"\eta~[m]", lab = L"η(t)")
plot!(tₘ⁺, ηₘ⁺, seriestype=:scatter, lab = "maxima")
display(plot!(tₘ¯, ηₘ¯, seriestype=:scatter, lab = "minima"))

# # Derivatives
# lbl_x = L"t~[s]";   lbl_y = "";   llbls = [L"\eta" L"\frac{d\eta}{dt}"]
# plt = plot(tₘ,[ηₘ η̇ₘ],xlab = lbl_x, ylab = lbl_y, lab = llbls)
# display(plt)
# llbls = [L"\eta" L"\frac{d\eta}{dt}"]
# plt = plot(tₘ,[η̇ₘ η̈ₘ], xlab = lbl_x, ylab = lbl_y, lab = llbls)
# display(plt)

# Spectra
display(nlnspec_plt)
display(spec_plt)
display(plot(scat[1], scat[2], xlabel = L"T_p~[s]", ylabel = L"H_s(T_p)~[m]", label = "Scatter diagram", linewidth=3))

#############################################################################################
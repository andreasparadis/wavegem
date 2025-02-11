# Load necessary modules
using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles
using FFTW, DSP, LsqFit, Statistics, BSplineKit
import ColorSchemes.darkrainbow

gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))
cb = darkrainbow

# Include necessary scripts for functions
include("text_process.jl")
include("signal_processing.jl")
include("peak_detect.jl")
include("wave_theory.jl")

###############################################################################
wout = false
# Select case
caseid = 4  # corresponds to the dname indice (0 for OW3D inp)
fileid = 2  # corresponds to the fname indice
dname = ("MaxFair", "MaxPitch", "MaxWave","Full")
fname = ("Nonlinear", "1st_order", "Gauss", "EWG_sum", "WG_fcsd")

###############################################################################
ProjPath::String = pwd()               # Project path
libpath = joinpath(ProjPath,"library") # Path to project's library folder

# Full scale simulation results
pdir::String = "SE" # Parent directory in library folder
case = joinpath("LCH","1")
fld = joinpath(libpath,pdir,case)
sub1 = joinpath(fld,"0")    # Subfolder 1
sub2 = joinpath(fld,"Decomposition")  # Subfolder 2
sub3 = joinpath(sub2,dname[caseid])   # Subfolder 3
sub4 = joinpath(sub3,"GR_Results")    # Subfolder 4

fid = ""

if caseid > 0 && caseid < 4 
    if fileid == 1
        fid = joinpath(sub3,"event") # Path to file
    elseif fileid == 2
        fid = joinpath(sub3,"event_lin") # Path to file
    elseif fileid == 3
        fid = joinpath(sub4,"etaG") # Path to file
    elseif fileid == 4
        fid = joinpath(sub4,"EWGsum") # Path to file
    elseif fileid == 5
        fid = joinpath(sub4,"etaFCSD") # Path to file
    end
elseif caseid == 4
    if fileid == 1
        fid = joinpath(sub2,"eta0") # Path to file
    elseif fileid == 2
        fid = joinpath(sub2,"eta_lin") # Path to file
    end
else
    fid = joinpath(sub1,"eta_t") # Path to file
end

cont = parse_fxw(fid, 1)

# Scaled results
λ = 100 # Scale factor (L_full/L_model)

# Read full scale results
fcut = 1/1.6 * sqrt(λ)

tᵢ = cont[:,1] ./ sqrt(λ) # time vector
ηᵢ = cont[:,2] ./ λ # Surface elevation
Lᵢ = length(tᵢ)
dtᵢ = (tᵢ[end]-tᵢ[1])/(Lᵢ-1)

# # Downsampling
# fₛₐₘₚ = 16
# dt = 1/fₛₐₘₚ
# L = Int(tᵢ[end]/dt + 1)
# t = range(0, tᵢ[end], L)
# itp_ηᵢ = interpolate(tᵢ, ηᵢ, BSplineOrder(4))
# η = itp_ηᵢ.(t)
# plot(tᵢ,ηᵢ)
# plot!(t,η)
###############################################################################
# Spectral analysis
# frᵢ, mgᵢ, phiᵢ, dfᵢ, Tₛᵢ, Nᵢfft, Hᵢ = one_side_asp(ηᵢ, tᵢ)
freq, mag, phi, df, Tₛ, Nfft, H = one_side_asp(ηᵢ, tᵢ)
fₚ = freq[findmax(mag)[2]];    ωₚ = 2π*fₚ
Tₚ = round(1/fₚ*1e3)/1e3
Hₛ = round(4*std(ηᵢ)*1e3)/1e3

# freq2 = fftfreq(Nfft,1/Tₛ)  # frequency range
# P2 = abs.(H)/Nfft        # two sided spectrum
# ϕ₂ = angle.(H)
###############################################################################
# Interpolation
Tᵣₑₜ = 128 # Return period
# Nf = max(512,Int(512/Tᵣₑₜ * round(tᵢ[end]))) # Number of frequency components
Nf = 512
fᶜ = 4.0    # New cut-off frequency
itp_sp = interpolate(freq, mag, BSplineOrder(4))
itp_phi = interpolate(freq, unwrap(phi), BSplineOrder(4))

fᵣₑₜ = 1/Tᵣₑₜ
# dF = (fᶜ-fᵣₑₜ)/(Nf-1)
# Nfi = Int(round(1/Tₛ/2/dF)) - 1
Fₛ = 1/((tᵢ[end]-tᵢ[1])/(Nf-1))

# FR = range(0,1/Tₛ/2, Nfi)
FR = range(fᵣₑₜ,fᵣₑₜ*Nf, Nf)
# FR = range(fᵣₑₜ,fᶜ, Nf)
MAG = itp_sp.(FR)
PHI = itp_phi.(FR);     PHI = unwrap(PHI)

dF = (FR[end]-FR[1])/(Nf-1)
τ = range(tᵢ[1],tᵢ[end],Nf)

# Reconstruct signal from interpolated spectrum
# FR2 = fftfreq(2*Nf-1,Fₛ)
H2 = Nf*MAG.*exp.(1im*PHI)
H2 = vcat(H2,conj(H2[end-1:-1:2]))

# ETA = zeros(ComplexF64,Nf)
# cₙ = MAG/2 .* exp.(1im*PHI)
# for i ∈ 1:Nf
#     for j ∈ 1:Nf
#         # ETA[j] = ETA[j] + H2[i]/Nf * exp(1im*2π*FR2[i]*τ[j])
#         # ETA[j] = ETA[j] + (cₙ[i]+conj(cₙ[i]))*exp(1im*2π*FR[i]*τ[j])
#         ETA[Nf+1-j] = ETA[Nf+1-j] + (cₙ[i]+conj(cₙ[i]))*exp(1im*2π*FR[i]*τ[j])*exp(1im*π)
#         # ETA[j] = ETA[j] + MAG[i] * cos(2π*FR[i]*τ[j]+PHI[i])
#     end
# end
# ETA = real(ETA)

# Hint = zeros(ComplexF64,Nfft)
# Hint[1:Nf] = Nf*MAG.*exp.(1im*PHI)
# Hint = Nf*MAG.*exp.(1im*PHI)

# Nf2 = Int(round(Nf/2))-1
# pad = zeros(Float64,Nf2)
# ETA = ifft([Hint; pad]) 
# # ETA = ifft(Hint)
# ETA = real(ETA[1:Nf])

# H2 = zeros(ComplexF64,Nfft)
# H2[1:Nf] = Hint 
# H2[Nfft-Nf+1:end] = conj(Hint[end:-1:1])
# H2 = vcat(Hint,conj(Hint[end-1:-1:2]))

ETA = real(ifft(H2))
τ = range(0,tᵢ[end],length(ETA))

###############################################################################
# Plots
plt_a = plot(xlab = L"t~[s]", ylab = L"\eta~[m]")
# plot!(tᵢ, ηᵢ, lab = L"H_s="*"$Hₛ [m]")
plot!(tᵢ, ηᵢ, lab="Original")
plot!(τ, ETA, lab="Reconstructed")

plt_b = plot(xlab = L"f~[Hz]", ylab = L"S(f)~[m^2 s]")
plot!(freq, mag, lab = L"T_p ="*"$Tₚ [s]")
plot!(FR,MAG, lab="Interpolated")
plot!(xlim=(0,fcut))

plt_c = plot(xlab = L"f~[Hz]", ylab = L"\phi~[rad]")
plot!(freq, unwrap(phi), lab="unwrap(ϕ)")
plot!(FR,PHI, lab="Interpolated")
plot!(xlim=(0,fcut))
# plot!(ylim=(-1.5*maximum(abs.(PHI[1:200])),1.5*maximum(abs.(PHI[1:200]))))

plt_span = plot(plt_a, plt_b, plt_c, layout = @layout [a ; c d])
display(plt_span)

###############################################################################
# Write output
if wout
    fld_out = joinpath(libpath,"UCL")
    sub_out = joinpath(fld_out,"Num_Targets")

    if caseid > 0 && caseid < 5
        ssub_out = joinpath(sub_out,dname[caseid])
        fid_out = joinpath(ssub_out,fname[fileid])
    else
        fid_out = joinpath(sub_out,"OW3D_inp")
    end

    open(fid_out, "w")
    writedlm(fid_out, [FR MAG PHI], '\t')

    if fileid == 2
        suffix = "HS00$(Int(1000*Hₛ))TP0$(Int(round(100*Tₚ)))TR$Tᵣₑₜ" # Suffix of output files
        fid_out = joinpath(ssub_out,dname[caseid]*".fronts")
        head0 = [suffix "" "" ""]
        head1 = ["f" "a" "angle" "phase"]
        head2 = ["Hz" "m" "rad" "rad"]
        open(fid_out, "w")
        writedlm(fid_out, [head0; head1; head2; FR MAG zeros(Float64,Nf) PHI], '\t')
    end

    if fileid == 4
        suffix = "HS00$(Int(1000*Hₛ))TP0$(Int(round(100*Tₚ)))TR$Tᵣₑₜ" # Suffix of output files
        fid_out = joinpath(ssub_out,dname[caseid]*"_EWGsum.fronts")
        head0 = [suffix "" "" ""]
        head1 = ["f" "a" "angle" "phase"]
        head2 = ["Hz" "m" "rad" "rad"]
        open(fid_out, "w")
        writedlm(fid_out, [head0; head1; head2; FR MAG zeros(Float64,Nf) PHI], '\t')
    end
end
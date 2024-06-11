using FFTW, DSP
using DelimitedFiles
using Plots, Statistics, LaTeXStrings, Dates
using BSplineKit

# Include necessary scripts for functions
include("signal_processing.jl")
include("jonswap.jl")
include("peak_detect.jl")
include("directories.jl")

#############################################################################################
# Read signal from text file in library folder
dir0::String = "SE/LAF/2/"
fid0::String = dir0*"Decomposition/"

#############################################################################################
# Import time-domain data for each pi/2 phase shifted realization
𝚶⃗ = Array{Float64}(undef,0,2)
# 1st column: Time    2nd column: Data
A00, A05, A10, A15 = (𝚶⃗[:] for _ = 1:4)
# A00
# fid ::String = "DS_Strath/A00_p.txt" # Path to file
fid ::String = fid0*"eta0" # Path to file
open("library/"*fid, "r")
A00 = readdlm("library/"*fid, '\t', Float64, '\n')
# A05
# fid ::String = "DS_Strath/A05_p.txt" # Path to file
fid ::String = fid0*"eta1" # Path to file
open("library/"*fid, "r")
A05 = readdlm("library/"*fid, '\t', Float64, '\n')
# A10
# fid ::String = "DS_Strath/A10_p.txt" # Path to file
fid ::String = fid0*"eta2" # Path to file
open("library/"*fid, "r")
A10 = readdlm("library/"*fid, '\t', Float64, '\n')
# A15
# fid ::String = "DS_Strath/A15_p.txt" # Path to file
fid ::String = fid0*"eta3" # Path to file
open("library/"*fid, "r")
A15 = readdlm("library/"*fid, '\t', Float64, '\n')

#############################################################################################
fcut = 0.625 # [Hz]

# ηᴴ = hilbert(A00[:,2]);   ηᴴₑₙᵥ = abs.(ηᴴ)
# A05 = A00[:,:];   A10 = A00[:,:];  A15 = A00[:,:] 
# A05[:,2] = real(ηᴴₑₙᵥ.*exp.(1im*(angle.(ηᴴ).+π/2)))
# A10[:,2] = real(ηᴴₑₙᵥ.*exp.(1im*(angle.(ηᴴ).+π)))
# A15[:,2] = real(ηᴴₑₙᵥ.*exp.(1im*(angle.(ηᴴ).+3π/2)))
# A05[:,2] .-= mean(A05);    A10[:,2] .-= mean(A10);    A15[:,2] .-= mean(A15)

# For signal truncation when needed
dt = A00[2,1]-A00[1,1]; tstart = 400; tend = tstart + 1600
ibeg = Int64(round(tstart/dt))+1;  iend = length(A00[:,1]) # Truncation limits
# ibeg = Int64(round(tstart/dt))+1;  iend = Int64(round(tend/dt))+1 # Truncation limits

#############################################################################################
# Refer Data values to 1st element of vector (make first Data element zero)
A̅00 = A00[ibeg:iend, 2]
A̅05 = A05[ibeg:iend, 2]
A̅10 = A10[ibeg:iend, 2]
A̅15 = A15[ibeg:iend, 2]

t00 = A00[ibeg:iend, 1] .- A00[ibeg, 1]
t05 = A05[ibeg:iend, 1] .- A05[ibeg, 1]
t10 = A10[ibeg:iend, 1] .- A10[ibeg, 1]
t15 = A15[ibeg:iend, 1] .- A15[ibeg, 1]

tOG = round.(t00*100)/100
nₜ = length(tOG)

# For JONSWAP spectrum for comparison at the end
Hₛ = 4*std(A̅00)
Tₑ = 25;     ω⁻ = 2π/Tₑ
Tᵢ = 1/fcut;    ω⁺ = 2π/Tᵢ;
dω = 2π/((nₜ-1)*dt)
Ncut = Int64(round((ω⁺-ω⁻)/dω))
#############################################################################################
# Obtain spectra for each realization via FFT
freq_1s, Â₀₀, ϕ₀₀, _, Ts, L = one_side_asp(A̅00, t00)
_, Â₀₅, ϕ₀₅, _, _, _ = one_side_asp(A̅05, t05)
_, Â₁₀, ϕ₁₀, _, _, _ = one_side_asp(A̅10, t10)
_, Â₁₅, ϕ₁₅, _, _, _ = one_side_asp(A̅15, t15)

# Spectral amplitudes (not scaled)
amp00 = Â₀₀*L
amp10 = Â₀₅*L
amp05 = Â₁₀*L
amp15 = Â₁₅*L

# Spectral phase angles
arg00 = unwrap(ϕ₀₀)
arg05 = unwrap(ϕ₀₅)
arg10 = unwrap(ϕ₁₀)
arg15 = unwrap(ϕ₁₅)

# Phase shifting
sp_peak = findmax(Â₀₀);     iₚ = sp_peak[2];    fₚ = freq_1s[iₚ];   Tₚ = 1/fₚ
fⱼ,Sⱼ = spectrum(Hₛ,Tₚ,3.3,Tᵢ,Tₑ,Ncut) # Generate JONSWAP spectrum for comparison at the end

N00 = phase_shift(arg00[iₚ])
arg00 .+= 2 * π * N00
N05 = phase_shift(arg05[iₚ])
arg05 .+= 2 * π * N05
N10 = phase_shift(arg10[iₚ])
arg10 .+= 2 * π * N10
N15 = phase_shift(arg15[iₚ])
arg15 .+= 2 * π * N15

# Decomposition
c00 = amp00 .* exp.(1im .* arg00)
c05 = amp05 .* exp.(1im .* arg05)
c10 = amp10 .* exp.(1im .* arg10)
c15 = amp15 .* exp.(1im .* arg15)

# Amplitudes and phase angles of each component
## Linear
cc1 = ((c00 - c10) - 1im * (c05 - c15)) / 4
AMP1, ARG1 = abs.(cc1), unwrap(angle.(cc1))
## 2nd+
cc2 = ((c00 + c10) - (c05 + c15)) / 4
AMP2, ARG2 = abs.(cc2), unwrap(angle.(cc2))
## 3rd
cc3 = ((c00 - c10) + 1im * (c05 - c15)) / 4
AMP3, ARG3 = abs.(cc3), unwrap(angle.(cc3))
## 2nd-
cc4 = (c00 + c05 + c10 + c15) / 4
AMP4, ARG4 = abs.(cc4), unwrap(angle.(cc4))

# Obtain time signals of the components via inverse FFT
L2 = Int(round(L/2))-1;  pad = zeros(Float64,L2)
t1 = ifft([cc1; pad]);  t1 = real(t1[1:nₜ])
t2 = ifft([cc2; pad]);  t2 = real(t2[1:nₜ])
t3 = ifft([cc3; pad]);  t3 = real(t3[1:nₜ])
t4 = ifft([cc4; pad]);  t4 = real(t4[1:nₜ])

total = real(t1+t2+t3+t4)

#############################################################################################
# Isolate, shift around 0 and plot 'pulses' based on peaks

## Find peaks and sort them in descending order
T_cut = 1/fcut

MinPeakVal = 2*std(A̅00)
MinPeakDist = T_cut
PeakVal, PeakPos, PeakId = peaks_max_ext(A̅00, tOG, MinPeakVal, MinPeakDist)

Val_LNR = real(t1)
MinPeakVal = 2*std(Val_LNR)
L_PeakVal, L_PeakPos, L_PeakId = peaks_max_ext(Val_LNR, tOG, MinPeakVal, MinPeakDist)

PId = sortperm(PeakVal, rev=true)
PeakVal = PeakVal[PId]
PeakPos = PeakPos[PId]
PeakId = PeakId[PId]

L_PId = sortperm(L_PeakVal, rev=true)
L_PeakVal = L_PeakVal[L_PId]
L_PeakPos = L_PeakPos[L_PId]
L_PeakId = L_PeakId[L_PId]

L_Peaks_len = length(L_PeakId)

## Isolate pulses given a specified time range (increase the range to consider more pusles)
t_range = 6  # [s]
rng = round(Int, (t_range/2) / Ts + 1)
NP = 8

for i ∈ 1:L_Peaks_len
    if L_PeakId[i] < rng+1
        L_PeakId[NP+1], L_PeakId[1] = L_PeakId[1], L_PeakId[NP+1]
        L_PeakPos[NP+1], L_PeakPos[1] = L_PeakPos[1], L_PeakPos[NP+1]
        L_PeakVal[NP+1], L_PeakVal[1] = L_PeakVal[1], L_PeakVal[NP+1]
    end
end

init0 = zeros(2*rng+1,NP)
t_pls, t_plsLNR, pulse, pulseLNR = (init0[:,:] for _ = 1:4)

## Isolate pulses around peaks
for i ∈ 1:NP
    # Time vector of t_range around peak i
    # Original signal
    t_pls[:, i] = (tOG[PeakId[i]-rng:PeakId[i]+rng]) .- tOG[PeakId[i]]
    pulse[:, i] = A̅00[PeakId[i]-rng:PeakId[i]+rng]
    
    # Linear signal
    t_plsLNR[:, i] = (tOG[L_PeakId[i]-rng:L_PeakId[i]+rng]) .- tOG[L_PeakId[i]]
    pulseLNR[:, i] = Val_LNR[L_PeakId[i]-rng:L_PeakId[i]+rng]
end

## Calculate mean pulse
timemean = t_pls[:,NP]
NLN_mean = mean!(ones(2*rng+1,1), pulse)

tmeanLNR = t_plsLNR[:,NP]
LNR_mean = mean!(ones(2*rng+1,1), pulseLNR)

#############################################################################################
# Compare pulses and calculate statistical measures
t_rng = 5  # [s]
id_rng = round(Int, (t_rng/2) / Ts + 1)  # +- in terms of time intervals
ID0 = div(size(pulseLNR, 1), 2)  # Index of t=0 (middle of vector)

## Pulses corresponding to the shortened time range
shrt = pulseLNR[ID0-id_rng:ID0+id_rng, :]

# Find the pulse with the highest peak
M_pls_i = findmax(shrt)[2]
fmf, inx = M_pls_i[1], M_pls_i[2]

## 'inx' is the identifier of the pulse containing the overall max
A = shrt[:, inx]

## Statistical comparison of other pulses against the identified pulse
## Exclude pulse A from the set of pulses
B = hcat(shrt[:, 1:inx-1], shrt[:, inx+1:end])

cvar = cov(A, B)
sigA = std(A)
sigB = std(B, dims=1)

## Pearson correlation coefficient
rho = cvar ./ (sigA * sigB)

R = A .- B
Ravg = mean(R, dims=1)
sigR = std(R, dims=1)

#############################################################################################
# Extract event
lb = PeakId[1]-500; ub = lb+1000 # For event plots

# itp = interpolate(tOG[lb:ub], A̅00[lb:ub], BSplineOrder(4))
# Lₑᵥ = nextpow(2,length(tOG[lb:ub])) * 2^7
# tₑᵥ = range(tOG[lb],tOG[ub], Lₑᵥ)
# Aₑᵥ = itp.(tₑᵥ)
# freq, mag, _,_,_,Nfft, H = one_side_asp(Aₑᵥ, tₑᵥ)
# FR, MAG, _,_,_,Nfft, H = one_side_asp(A̅00[lb:ub], tOG[lb:ub])

Aₑᵥ = zeros(Float64, nₜ)
Aₑᵥ[lb:ub] = A̅00[lb:ub]
freq, mag, _ = one_side_asp(Aₑᵥ, tOG)
FR, MAG, _ = one_side_asp(A̅00[lb:ub], tOG[lb:ub])

# itp_sp = interpolate(freq, mag, BSplineOrder(2))
# FR = range(0,fcut, Nfft)
# MAG = itp_sp.(FR)
# EV = ifft(H)

plt = plot(freq, mag, yscale=:log10)
plot!(FR, MAG, yscale=:log10)
plot!(xlim=(0,fcut))
display(plt)
#############################################################################################
# Surface elevation (temporal)
fid_elev::String = "eta_lin" # File name
open("library/"*fid0*fid_elev, "w")
writedlm("library/"*fid0*fid_elev, [tOG t1], '\t')

# Surface elevation spectrum
fid_elev::String = "eta_lin_spec" # File name
open("library/"*fid0*fid_elev, "w")
writedlm("library/"*fid0*fid_elev, [freq_1s AMP1], '\t')

# Event
evdir::String = "events/"
# if isdir(evdir)
#     error("Andreas: This folder already exists, since this event has already been recorded!")
# else
#     mkdir("library/"*fid0*evdir)
# end

fid_ev::String = "event_1" # File name
open("library/"*fid0*evdir*fid_ev, "w")
writedlm("library/"*fid0*evdir*fid_ev, [tOG[lb:ub] A̅00[lb:ub]], '\t')

# lbe = Int(round(360/dt)); ube = Int(round(415/dt))
# fid_ev::String = "event_2" # File name
# open("library/"*fid0*evdir*fid_ev, "w")
# writedlm("library/"*fid0*evdir*fid_ev, [tOG[lbe:ube] A̅00[lbe:ube]], '\t')

fid_ev::String = "event_1_int" # File name
open("library/"*fid0*evdir*fid_ev, "w")
writedlm("library/"*fid0*evdir*fid_ev, [tOG Aₑᵥ], '\t')

#############################################################################################
## PLOTS 
fid_fig::String = "library/"*fid0*"/Figures/"
if !isdir(fid_fig)
    mkdir("library/"*fid0*"/Figures")
end

# Plot decomposed signal and its components
plt1 = plot(tOG[lb:ub], t1[lb:ub], xlab = "Time [sec]", ylab = "Elevation [m]", line = :solid, color = :blue, label = "Linear")
plot!(tOG[lb:ub], t2[lb:ub], line = :solid, color = :green, label = "2nd+")
plot!(tOG[lb:ub], t3[lb:ub], line = :solid, color = :magenta, label = "3rd")
plot!(tOG[lb:ub], t4[lb:ub], line = :solid, color = :black, label = "2nd-")
plot!(tOG[lb:ub], A̅00[lb:ub], line = :dash, color = :red, linewidth = 1, label = "Original")
display(plt1)
savefig(fid_fig*"decomp_elev.svg")
savefig(fid_fig*"decomp_elev.png")

plt2 = plot(tOG[lb:ub], [total[lb:ub] A̅00[lb:ub]], xlab = "Time [sec]", ylab = "Elevation [m]", line = [:solid :dash], color = [:blue :red], label = ["Sum" "Original"])
savefig(fid_fig*"elev_compare.svg")
savefig(fid_fig*"elev_compare.png")

# plt3 = plot(tₑᵥ, Aₑᵥ, line = :solid, color = :blue, linewidth = 1, label = "Interpolation")
# plot!(tOG[lb:ub], A̅00[lb:ub], line = :dash, color = :red, linewidth = 1, label = "Original")
# display(plt3)

# Plot identified peaks on top original and linear signals
plot(tOG, A̅00, xlab = "Time [sec]", ylab = "Elevation [m]",  label = "Original signal")
display(plot!(PeakPos, PeakVal, seriestype=:scatter, ms=2, mc=:red, lab = "Peaks",legend=:outerbottom, legendcolumns=2))
savefig(fid_fig*"elev_peaks_nln.svg")
savefig(fid_fig*"elev_peaks_nln.png")

plot(tOG, Val_LNR, xlab = "Time [sec]", ylab = "Elevation [m]",  label = "Linear signal")
display(plot!(L_PeakPos, L_PeakVal, seriestype=:scatter, ms=2, mc=:red, lab = "Peaks",legend=:outerbottom, legendcolumns=2))
savefig(fid_fig*"elev_peaks_lin.svg")
savefig(fid_fig*"elev_peaks_lin.png")

# Plot 'pulses' and mean
plot(t_pls, pulse, title="Nonlinear Signal", xlabel="Time [sec]", ylabel="Elevation [m]", grid=true)
plta = plot!(timemean,NLN_mean,line = :dash, color = :red, linewidth = 2, label = "Mean")

plot(t_plsLNR, pulseLNR, title="Linearized Signal", xlabel="Time [sec]", ylabel="Elevation [m]", grid=true)
pltb = plot!(tmeanLNR,LNR_mean,line = :dash, color = :red, linewidth = 2, label = "Mean")

display(plot(plta, pltb, layout = @layout [a ; b ]))
savefig(fid_fig*"elev_pulses_compare.svg")
savefig(fid_fig*"elev_pulses_compare.png")

# Plot spectral decomposition
df = (freq_1s[end]-freq_1s[1])/(L/2);
dfⱼ = abs(fⱼ[end]-fⱼ[1])/(length(fⱼ)-1)
η̂ = sqrt.(2*Sⱼ*dfⱼ)

flb = 1; fub = Int(round(fcut/df))
plt_sp = plot(freq_1s[flb:fub], Â₀₀[flb:fub], xlab = "f [Hz]", ylab = L"S(f)~[m^2 s]", line =:dashdot, lab="Original", yscale=:log10)
plot!(freq_1s[flb:fub], AMP1[flb:fub]/L, line =:solid, lw=0.5, label = "Linear", yscale=:log10)
plot!(freq_1s[flb:fub], AMP2[flb:fub]/L, line =:dash, lw=0.5, label = "2nd+", yscale=:log10)
plot!(freq_1s[flb:fub], AMP3[flb:fub]/L, line =:dot, lw=0.5, label = "3rd", yscale=:log10)
plot!(freq_1s[flb:fub], AMP4[flb:fub]/L, line =:solid, lw=0.5, label = "2nd-", yscale=:log10)
plot!(fⱼ, η̂, lab="JONSWAP", color=:red, lw=2, yscale=:log10)
plot!(xlim=(0,fcut), ylim=(1e-6,1), minorgrid=true)
display(plt_sp)
savefig(fid_fig*"decomp_spectra.svg")
savefig(fid_fig*"decomp_spectra.png")

freq, mag, _ = one_side_asp(Aₑᵥ,tOG)
plt_ev = plot(freq_1s[flb:fub], Â₀₀[flb:fub], xlab = "f [Hz]", ylab = L"S(f)~[m^2 s]", line =:dashdot, lab="Original", yscale=:log10)
plot!(freq, mag, line =:solid, lab="Event", yscale=:log10)
plot!(fⱼ, η̂, lab="JONSWAP", color=:red, lw=2, yscale=:log10)
plot!(xlim=(0,fcut), ylim=(1e-6,1), minorgrid=true)
display(plt_ev)
savefig("library/"*fid0*evdir*"event_spectra.svg")
savefig("library/"*fid0*evdir*"event_spectra.png")

# Compare spectra of input signal and OW3D simulation
elev_fid::String = "library/"*dir0*"0/eta_t"
open(elev_fid, "r")   
elev = readdlm(elev_fid, '\t', Float64, skipstart=1, '\n')
tinp, ηinp = elev[:,1], elev[:,2]
freq, mag, _ = one_side_asp(ηinp,tinp)

plt_comp = plot(freq, mag, xlab = "f [Hz]", ylab = L"S(f)~[m^2 s]", lab="Input", yscale=:log10)
plot!(freq_1s[flb:fub], Â₀₀[flb:fub], line =:dot, lab="Original", yscale=:log10)
plot!(freq_1s[flb:fub], AMP1[flb:fub]/L, line =:dashdot, lw=0.2, label = "Linear", yscale=:log10)
plot!(fⱼ, η̂, lab="JONSWAP", color=:red, lw=2, yscale=:log10)
plot!(xlim=(0,fcut), ylim=(1e-6,1), minorgrid=true)
display(plt_comp)
savefig(fid_fig*"spect_comp.svg")
savefig(fid_fig*"spect_comp.png")

h = PeakVal[1]
plt_groups = plot(tOG, A̅00, xlab = "Time [sec]", ylab = "Elevation [m]",  label = "Original signal")
for i ∈ 1:max(20, length(PeakVal))
    c = PeakPos[i]
    w = 25
    plot!(Shape([c-w,c-w,c+w,c+w],[-h,h,h,-h]), color=:red, opacity=0.2, label = "WG$i")
end
plot!(legend = false)
display(plt_groups)
savefig(fid_fig*"WGs.svg")

plt_groups = plot(tOG, A̅00, xlab = "Time [sec]", ylab = "Elevation [m]",  label = "Original signal")

c = PeakPos[1]
w = 50
plot!(Shape([c-w,c-w,c+w,c+w],[-h,h,h,-h]), color=:red, opacity=0.3, label = "Event")
plot!(xlim=(6000,7000))
display(plt_groups)
savefig(fid_fig*"Event1.svg")

# # Find zero Up or Down Crossings
# UCiD = findall(diff(sign.(A̅00)) .== 2)  # Up-crossing (-2 for Down-crossing)
# UC_Val = A̅00[UCiD]
# UC_Pos = tOG[UCiD]

# UCiD = findall(diff(sign.(Val_LNR)) .== 2)  # Up-crossing (-2 for Down-crossing)
# L_UC_Val = Val_LNR[UCiD]
# L_UC_Pos = tOG[UCiD]

# # Plot original signal along with up crossing points
# plt = plot(tOG, A̅00, line=:solid, color=:blue, lab="Elevation")
# plot!(UC_Pos, UC_Val, seriestype=:scatter, mc=:red, ms=2, lab = "Up-Crossing points")
# plot!(
#     xlim=(tstart, tstart+200),
#     xlabel="Time [sec]", ylabel="Elevation [m]",
#     legend=:topright
# )
# display(plt)
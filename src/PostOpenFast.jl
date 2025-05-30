module PostOpenFast
using Plots, LaTeXStrings, DelimitedFiles, Dates, LinearAlgebra
using FFTW, DSP, Statistics, BSplineKit, StatsBase, Distributions

include("WAVEGEM.jl")
import .WAVEGEM

# Include necessary scripts for functions
include("func/signal_processing.jl")
include("func/text_process.jl")
include("func/peak_detect.jl")
include("func/jonswap.jl")
include("func/event_proc.jl")

import ColorSchemes.darkrainbow
gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))
cb = darkrainbow
#############################################################################################
# Import Global Variables & Module specific inputs
ρ, g, _, d, _, _, γ, fcut = WAVEGEM.GlobInp0
runstr, case_id, casedir, rundir, phipath, OW3Dcdir, OW3Drdir, OW3Dphipath, 
            Decpath, DecEvs, DecFigs, postOFpath, postOW3Dpath, OFASTpath = WAVEGEM.GlobPaths
Wave = WAVEGEM.Wave
full, fev = WAVEGEM.POFflags
CET, evID = WAVEGEM.CET, WAVEGEM.evID
T₀ₛ, T₀ₕ, T₀ₚ = WAVEGEM.FOWT[:]

#############################################################################################
# Open and read files
## Full simulation
OFASTout = joinpath(postOFpath,"outD_")   # OpenFAST output file
cont = parse_fxw(OFASTout, 5)
heads = readdlm(joinpath(postOFpath,"dataHdr"))

## Surface elevations for comparative plots
sefname = joinpath(Decpath,"eta_nln")      # Surface elevation from OW3D (0 phase shift)
ηOW3D = parse_fxw(sefname, 0)

#############################################################################################
# Assign Variables
NoRows = size(cont)[2]

# Nₜ = size(cont)[1]
# dt = 0.01
# t = zeros(Float64,Nₜ)
# [t[i] = (i-1)*dt for i ∈ 1:Nₜ]  

# Surge, Heave, Pitch, Wave1Elev, Wave1Elev1, FAIRTEN1, FAIRTEN2, AnchTen1, AnchTen2,
# HydroFxi, HydroFzi = [cont[:,i] for i = 1:NoRows]

# t, Surge, Heave, Pitch, Wave1Elev, Wave1Elev1, FAIRTEN1, FAIRTEN2, 
# HydroFxi, HydroFzi, HydroMyi = [cont[:,i] for i = 1:NoRows]

t, Surge, Heave, Pitch, Wave1Elev, Wave1Elev1, FAIRTEN1, FAIRTEN2, 
HydroFxi, HydroFzi, HydroMyi, PtfmTAxt, PtfmTAzt = [cont[:,i] for i = 1:NoRows]

aCoM = sqrt.(PtfmTAxt.^2 .+ PtfmTAzt.^2)  # CoM acceleration
Nₜ = length(t)       
dt = t[2]-t[1]

#############################################################################################
# Critical Event statistics & selection
RespVar = zeros(Float64,Nₜ) # The response variable to be analysed

if CET == 1        # Fairlead tension critical event
    RespVar[:] = FAIRTEN2[:]
    case_str = "MaxFair"
    varname = "Fairlead tension"
elseif CET == 2    # Pitch critical event
    RespVar[:] = Pitch[:]
    case_str = "MaxPitch"
    varname = "Pitch"
elseif CET == 3    # CoM extreme displacement event 
    RespVar[:] = disp[:]
    case_str = "MaxCoM"
    varname = "Centre of mass"
else                # Response to extreme wave event
    RespVar[:] = Wave1Elev[:]
    case_str = "MaxWave"
    varname = "Surface elevation"
end
figOF = joinpath(postOFpath,case_str)
## Make directory for specified event if it doesn't exist
if !isdir(figOF)
    mkdir(figOF)
end

#############################################################################################
# Response Statistics
StatMetrics = [mean(RespVar); std(RespVar); median(RespVar); rms(RespVar); skewness(RespVar); kurtosis(RespVar)-3]
RefVar = StatMetrics[4]   # Reference/normalisation variable (rms, median, mean, RespVar[1])
## Find all maxima above a given threshold of the selected response time history
RespVar_Peaks, _ = peaks_max_ext(RespVar.-RefVar,t, 0, 0, false)

## Distribution of normalised peaks with reference to the rms value - Weibull
SampVar = RespVar_Peaks./RefVar   # Sample
WeibFit = fit(Weibull, SampVar);  kʷ = WeibFit.α;    λʷ = WeibFit.θ # Distribution fit
WeibRange = range(minimum(SampVar), maximum(SampVar), 2^9);  WeibPDF = pdf.(WeibFit, WeibRange)  # Theoretical distribution
WeibCDF = cdf.(WeibFit, WeibRange)
ptile = 1/100;  F̅₁ₚ = λʷ * (-log(ptile))^1/kʷ;     F₁ₚ = (F̅₁ₚ+1)*RefVar
### Plot PDF - Compare sample & fit
hplt_F2 = histogram(SampVar, normalize=:pdf, lab=:false, palette=[cb[8]; cb[11]; cb[2]; cb[5]])
plot!(hplt_F2, title="PDF - $varname maxima", xlab=L"\frac{F_p}{F^{rms}}-1", ylab="P(f)")
plot!(hplt_F2, WeibRange, WeibPDF, lw=2, lab="Weibull: κ=$(round(kʷ*1e3)/1e3), λ=$(round(λʷ*1e3)/1e3)")
plot!(hplt_F2,[mean(SampVar)+3*std(SampVar); mean(SampVar)+3*std(SampVar)], [0; maximum(WeibPDF)], ls=:dashdot, lab=L"\mu+3\sigma")
plot!([F̅₁ₚ; F̅₁ₚ], [0; maximum(WeibPDF)], ls=:dashdot, lab="Q($(ptile*100)%)")
plot!(twinx(),WeibRange, 1 .-WeibCDF, lw=2, ylab="Quantile", lab="Q(p)", ylim=(0,1), legendposition=:bottomright, palette=[cb[5]], ls=:dash)
### Q-Q plot of sample & fit
data_quantiles = sort(SampVar)
theor_quantiles = quantile.(WeibFit, (1:length(SampVar)) / (length(SampVar) + 1))
QQ_WeibFit = plot(xlabel="Theoretical Quantiles", ylabel="Sample Quantiles", title="Q-Q Plot - $varname maxima", palette=[cb[11]; cb[8]])
plot!(QQ_WeibFit, theor_quantiles, data_quantiles, lab="Weibull: κ=$(round(kʷ*1e3)/1e3), λ=$(round(λʷ*1e3)/1e3)", lw=2)
plot!(QQ_WeibFit, theor_quantiles, theor_quantiles, linestyle=:dash, label="45°")

## Distribution of signal values normalised to the rms value (Log-Normal)
SampVar = RespVar./RefVar # Sample
LogNormFit = fit(LogNormal, SampVar);    μˡⁿ = LogNormFit.μ;    σˡⁿ = LogNormFit.σ   # Distribution fit
LogNormRange = range(0, maximum(SampVar), 2^9); LogNormPDF = pdf.(LogNormFit,LogNormRange)    # Theoretical distribution
### Plot PDF - Compare sample & fit
hplt_LNF2 = histogram(SampVar, normalize=:pdf, lab=:false, palette=[cb[8]; cb[11]; cb[5]])
plot!(hplt_LNF2, title="PDF - $varname response", xlab=L"\frac{F}{F_{rms}}", ylab="P(f)")
plot!(hplt_LNF2, LogNormRange, LogNormPDF, lw=2, lab="Log-Normal: μ=$(round(μˡⁿ*1e3)/1e3), σ=$(round(σˡⁿ*1e3)/1e3)")
## Q-Q plot
data_quantiles = sort(SampVar)
theor_quantiles = quantile.(LogNormFit, (1:length(SampVar)) / (length(SampVar) + 1))
QQ_LogNorm = plot(xlabel="Theoretical Quantiles", ylabel="Sample Quantiles", title="Q-Q Plot - $varname", palette=[cb[11]; cb[8]])
plot!(QQ_LogNorm, theor_quantiles, data_quantiles, lab="μ=$(round(μˡⁿ*1e3)/1e3), σ=$(round(σˡⁿ*1e3)/1e3)", lw=2)
plot!(QQ_LogNorm, theor_quantiles, theor_quantiles, linestyle=:dash, lab="45°")

StatMetrics = round.(StatMetrics.*1e3)./1e3
fstmet = joinpath(postOFpath,case_str,"StatMetrics")
Chead = ["μ" "σ" "med()" "rms()" "skew()" "kurt()"]
open(fstmet, "w")
writedlm(fstmet, [Chead;StatMetrics'], '\t')

## Find and sort extreme maxima in the selected response time history
# cstd = 5     # Standard deviation coeffcient  (default=5)
# MinPeakVal = cstd*std(RespVar) .+ mean(RespVar)    # Amplitude threshold
MinPeakVal = F₁ₚ    # Amplitude threshold
cstd = round((MinPeakVal-StatMetrics[1])/StatMetrics[2]*1e3)/1e3     # Standard deviation coeffcient
Aᵢ, tᵢ, i⁺, _ = peaks_max_ext(RespVar,t, MinPeakVal, 0, true)
Nce = length(Aᵢ)        # Number of peaks above threshold

# Identify corresponding wave events
ωᴴ,_ = instant_freq(Wave1Elev,t)    # Instantaneous frequency
## Starting point of corresponding event
lbs = zeros(Int64,Nce)
iᶠ = Array{Int64}(undef,0)
Aᶠ = Array{Float64}(undef,0)
tᶠ = Array{Float64}(undef,0)

for i ∈ 1:Nce
    ilb = i⁺[i]
    while sign(ωᴴ[ilb]) == sign(ωᴴ[ilb-1]) && ilb > 2
        ilb = ilb - 1
    end
    # Remove duplicates
    if !(ilb ∈ lbs)
        push!(iᶠ,i⁺[i])
        push!(Aᶠ, Aᵢ[i])
        push!(tᶠ, tᵢ[i])
    end
    lbs[i] = ilb
end

### Store time instances of maxima
f_tinst = joinpath(postOFpath,case_str*"_tinsts")
if cstd != 5
    f_tinst = f_tinst*"_$(Int(ptile*100))ptile"
end
open(f_tinst, "w")
writedlm(f_tinst, [tᶠ Aᶠ], '\t')

# Statistical metrics for all responses
RespStats = zeros(Float64,NoRows-1,7)   # Matrix summary of response statistics
for i ∈ 2:NoRows
    RVar = cont[:,i]
    RespStats[i-1,1] = mean(RVar)
    RespStats[i-1,2] = std(RVar)
    RespStats[i-1,3] = median(RVar)
    RespStats[i-1,4] = rms(RVar)
    RespStats[i-1,5] = skewness(RVar)
    RespStats[i-1,6] = kurtosis(RVar)-3
    RespStats[i-1,7] = maximum(RVar)
end
RespStats = round.(RespStats.*1e3)./1e3
fstats = joinpath(postOFpath,"RespStats")
Chead = ["μ" "σ" "med()" "rms()" "skew()" "kurt()" "max()"]
Rhead = ["Response"; heads[2:NoRows]]
open(fstats, "w")
writedlm(fstats, [Rhead [Chead;RespStats]], '\t')

# Time histories
if fev
    # Event
    cont_ev = parse_fxw(figOF*"/outD_EV$evID", 5)
    # FWG
    cont_FWG = parse_fxw(figOF*"/outD_EV$(evID)_ReFoGWs", 5)
    # # DAM
    # cont_DAM = parse_fxw_pf(figOF*"/outD_DAM", 0, 5)
    # # 2AM
    # cont_2AM = parse_fxw_pf(figOF*"/outD_2AM", 0, 5)
    # # ALT_2AM
    # cont_ALT_2AM = parse_fxw_pf(figOF*"/outD_ALT_2AM", 0, 5)
    # # SFWG
    # # cont_SFWG = parse_fxw_pf(figOF*"/outD_SFWG", 0, 5)
    # # DWG
    # cont_DWG = parse_fxw_pf(figOF*"/outD_DWG", 0, 5)
    # # BEAT
    # cont_BEAT = parse_fxw_pf(figOF*"/outD_BEAT", 0, 5)
    # # NewWave
    # cont_NW = parse_fxw_pf(postOFpath*"/outD_NW", 0, 5)
    trange = Int(round(cont_ev[end,1]/2/dt))
    imax = findmax(cont_ev[:,8])[2] # Position of max fairten  
else
    trange = Int(round(2^7/dt))
    imax = trange
end

# Truncated time vector for event plots 
peak = [tᶠ[evID];iᶠ[evID]]    # Time and interval of event peak
lb = Int(peak[2]) - imax
if lb < 1
    lb=1
end
ub = Int(lb+2*trange-1) 
tₑᵥ = t[lb:ub];     Lₑᵥ = length(tₑᵥ)

#############################################################################################
# PLOTS
## Plot all stored variables of the whole simulation
if full
    # Plot OpenFast results
    for i ∈ 2:NoRows
        plti = plot(xlab = "t [s]", title = "$(heads[i])", legend=:topleft, palette=[cb[11]])
        plot!(t, cont[:,i], lw=2, lab = "Full sim")
        display(plti)
        savefig(joinpath(postOFpath,"$(heads[i]).svg"))
        savefig(joinpath(postOFpath,"$(heads[i]).png"))
    end

    # Wave Spectra
    FR, MAG, _ = one_side_asp(Wave1Elev, t)               # OpenFAST - spectrum of surface elevation
    FR_1st, MAG_1st, _ = one_side_asp(Wave1Elev1, t)      # OpenFAST - spectrum of 1st order elevation
    FR_2nd, MAG_2nd, _ = one_side_asp(Wave1Elev .- Wave1Elev1, t)   # OpenFAST - spectrum 2nd order
    FR_OW3D, MAG_OW3D, _ = one_side_asp(ηOW3D[:,2], ηOW3D[:,1])     # Spectrum of OW3D elevation (nln)
    fₚ = FR[findmax(MAG)[2]]

    # Plot OpenFAST surface elevation spectra
    plt = plot(xlab = L"f~[Hz]", ylab = L"S(f)~[m]", palette=[cb[11];cb[4];cb[8]]) 
    plot!(FR, MAG, yscale=:log10, line =:dot, lab = "OpenFAST")
    plot!(FR_1st, MAG_1st, yscale=:log10, opacity=0.8, lab = L"1^{st}~\mathrm{order}")
    plot!(FR_2nd, MAG_2nd, yscale=:log10, opacity=0.6, lab = L"2^{nd}~\mathrm{order}")
    plot!(xlim=(0,fcut), ylim=(1e-6,1),  minorgrid=true)
    display(plt)
    savefig(joinpath(postOFpath,"OFWaveSpectra.svg"))
    savefig(joinpath(postOFpath,"OFWaveSpectra.png"))

    # Comparative plots of surface elevation spectra: OpenFAST vs OW3D 
    plt = plot(xlab = L"f~[Hz]", ylab = L"S(f)~[m]", legendcolumns=:2, palette=[cb[11];cb[4];cb[8]]) 
    plot!(FR_OW3D, MAG_OW3D, lab = "OW3D")
    plot!(FR, MAG, line =:dash, lw = 0.5, lab = "OpenFAST")
    plot!(1/T₀ₛ*ones(Float64,100), maximum(MAG)*range(0,1,100), line=:dash,  lw = 2, lab = L"f_{0}^{surge}")
    plot!(1/T₀ₕ*ones(Float64,100), maximum(MAG)*range(0,1,100), line=:dash, lw = 2, lab = L"f_{0}^{heave}")
    plot!(1/T₀ₚ*ones(Float64,100), maximum(MAG)*range(0,1,100), line=:dash, lw = 2, lab = L"f_{0}^{heave}")
    plot!(xlim=(0, 0.625))
    display(plt)
    savefig(joinpath(postOFpath,"CompWaveSpectra.svg"))
    savefig(joinpath(postOFpath,"CompWaveSpectra.png"))
end

## Comparative multi-plot of key responses to selected extreme event
### Subplot 1 - Surge and CoM displacement
plt_eva = plot(legend=:topleft, palette=[cb[10];cb[10]])
plot!(tₑᵥ, Surge[lb:ub], lw=2, ylab="[m]", lab="Surge")
# plot!(tₑᵥ, aCoM[lb:ub], lw=1, lab=L"a_{CM}", line=:dot)
plot!(tₑᵥ, maximum(Surge)*ones(Float64,Lₑᵥ), lab="max", line=:dash)
# plot!(tₑᵥ, maximum(aCoM)*ones(Float64,Lₑᵥ), lab="max", line=:dot, lw=1)
### Subplot 2 - Pitch angle
plt_evb = plot(legend=:topleft, palette=[cb[4];cb[4]])
plot!(tₑᵥ, Pitch[lb:ub], lw=2, ylab="[deg]", lab="Pitch")
plot!(tₑᵥ, maximum(Pitch)*ones(Float64,Lₑᵥ), lab="max", line=:dash)
### Subplot 3 - Fairlead tension
plt_evc = plot(xlab="t [s]", legend=:topleft, palette=[cb[9];cb[9]])
plot!(tₑᵥ, FAIRTEN2[lb:ub], lw=2, ylab="[kN]", lab = "Fairlead")
plot!(tₑᵥ, maximum(FAIRTEN2)*ones(Float64,Lₑᵥ), lab="max", line=:dash)
### Subplot 4 - Surface elevation and Heave response
plt_evd = plot(xlab="t [s]", legend=:bottomleft, legendcolumns=:2, palette=[cb[1];cb[8];cb[1];cb[8]])
plot!(tₑᵥ, Wave1Elev[lb:ub], lw=2,  lab=L"\eta")
# plot!(tₑᵥ, Wave1Elev1[lb:ub], lw=2,  lab=L"\eta^1")
plot!(tₑᵥ, Heave[lb:ub], lw=2, ylab="[m]", lab="Heave")
plot!(tₑᵥ, maximum(Wave1Elev)*ones(Float64,Lₑᵥ), lab=L"max(\eta)", line=:dash)
# plot!(tₑᵥ, maximum(Wave1Elev1)*ones(Float64,Lₑᵥ), lab=L"max(\eta^1)", line=:dash)
plot!(tₑᵥ, maximum(Heave)*ones(Float64,Lₑᵥ), lab=L"max(z)", line=:dash)
### Concentrated plot (1+2+3+4)
plt_event = plot(plt_eva, plt_evb, plt_evc, plt_evd, layout = @layout [a b; c d])
display(plt_event)
savefig(joinpath(figOF,"EE$(evID)_resp.svg"))
savefig(joinpath(figOF,"EE$(evID)_resp.png"))

# Plot the selected response variable
pltRV = plot(t, RespVar, xlab="t [s]", ylab="RespVar", title=case_str, lab=:false, palette=[cb[11];cb[4];cb[8]])
plot!([t[1],t[end]], [MinPeakVal,MinPeakVal], line=:dot, lw=2, lab=L"\mu+"*"$(cstd)std")
plot!([0;t[end]],[rms(RespVar);rms(RespVar)], ls=:dash, lw=2, lab="rms(RespVar)")
display(pltRV)
# Distribution of RespVal peaks
## Weibull fit & Log-normal fit
plt_DistF2_1 = plot(hplt_F2, QQ_WeibFit, layout=@layout [a;b])
display(plt_DistF2_1)
savefig(joinpath(figOF,"Dist_FT2_max.svg"))
savefig(joinpath(figOF,"Dist_FT2_max.png"))
plt_DistF2_2 = plot(hplt_LNF2, QQ_LogNorm, layout=@layout [a;b])
display(plt_DistF2_2)
savefig(joinpath(figOF,"Dist_FT2.svg"))
savefig(joinpath(figOF,"Dist_FT2.png"))

#############################################################################################
if fev
    for i ∈ 2:NoRows
        plti = plot(xlab = "t [s]", title = "$(heads[i])", legend=:topleft, palette=[cb[11];cb[4];cb[8];cb[1]])
        plot!(tₑᵥ, cont[lb:ub,i], lw=3, lab = "Full sim")
        plot!(tₑᵥ, cont_ev[:,i], lw=2, lab = "Event")
        plot!(tₑᵥ, cont_FWG[:,i], line=:dot, lw=2, lab = L"ReFoGWs")
        # plot!(tₑᵥ, cont_SFWG[:,i], line=:dot, lw=2, lab = L"\sum g_n (t)")
        # plot!(tₑᵥ, cont_NW[:,i], line=:dot, lw=2, lab = "Focused Wave")
        display(plti)
        # savefig(figOF*"/$(heads[i]).svg")
        # savefig(figOF*"/$(heads[i]).png")
    end

    pltFairA = plot(xlab = "t [s]", title = "$(heads[8])", legend=:topleft, palette=[cb[11];cb[4];cb[8];cb[1]])
    plot!(tₑᵥ, FAIRTEN2[lb:ub], lw=3, lab = "Full sim")
    plot!(tₑᵥ, cont_FWG[:,8], lw=2, lab = L"ReFoGWs")
    display(pltFairA)
    savefig(joinpath(figOF,"FairA_ev$evID.svg"));  savefig(joinpath(figOF,"FairA_ev$evID.png"))
end
end

# Move to response statistics to examine further responses (example here is Surge response)
# ## Surge
# μS = mean(Surge);    σS = std(Surge)
# Nrang = range(μS-4*σS, μS+4*σS, 2^9);  Ncurve = 1/(σS*sqrt(2π)) * exp.(-((Nrang.-μS)/(sqrt(2)*σS)).^2)
# hplt_S = histogram(Surge, normalize=:pdf, lab=:false)
# plot!(hplt_S, title="Probability density function of Surge", xlab=L"Surge~[m]", ylab=L"P(Surge)")
# plot!(hplt_S, Nrang, Ncurve, lw=2, lab="Normal Fit: μ=$(round(μS*1e3)/1e3), σ=$(round(σS*1e3)/1e3)")
# ### Q-Q plot
# NormS = fit(Normal,Surge)
# data_quantiles = sort(Surge)
# theor_quantiles = quantile.(NormS, (1:length(Surge)) / (length(Surge) + 1))
# plot(theor_quantiles, data_quantiles, xlabel="Theoretical Quantiles", ylabel="Sample Quantiles", title="Surge Q-Q Plot", lab="μ=$(round(NormS.μ*1e3)/1e3), σ=$(round(NormS.σ*1e3)/1e3)")
# plot!(theor_quantiles, theor_quantiles, linestyle=:dash, color=:red, label="45-degree line")

# # Surge low-pass filter
# fcˡᵖ = fcut/32;     fₛˡᵖ = round(1/dt);     Ntaps = nextpow(2,fₛˡᵖ/fcˡᵖ)   # No of taps
# Surge_Low = low_pass_filter(Surge,fcˡᵖ,fₛˡᵖ,Ntaps)
# display(plot(t, Surge_Low, title="Surge-Low-pass filter"))
# # frS, magS,_ = one_side_asp(Surge_Low.-mean(Surge_Low),t)
# # plot(frS,magS)
# # plot!(xlim=(0,4*fcˡᵖ))
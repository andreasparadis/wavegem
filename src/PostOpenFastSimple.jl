module PostOpenFastSimple
using Plots, LaTeXStrings, DelimitedFiles, Dates, LinearAlgebra
using FFTW, DSP, Statistics, BSplineKit, StatsBase, Distributions, StatsPlots

include("WAVEGEM.jl")
import .WAVEGEM

# Include necessary scripts for functions
include("func/signal_processing.jl")
include("func/text_process.jl")
include("func/peak_detect.jl")
include("func/jonswap.jl")
include("func/event_proc.jl")
include("func/stat_process.jl")

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
res_path = joinpath(postOFpath,"MaxFair")
res_dir = readdir(res_path)
HG_fcount = count(file -> startswith(file, "outD_HG"), res_dir)
heads = readdlm(joinpath(res_path,"dataHdr"))
Øᵃⁿʸ = Array{Any}(undef,0)
AllPlots, plt_img_names = [Øᵃⁿʸ[:] for _ = 1:2]
∅ᴴᴳ = zeros(Float64,HG_fcount)
Fₘₐₓ, tₘₐₓ = [∅ᴴᴳ[:] for _ = 1:2]

for i ∈ 1:HG_fcount
    # Open and read files
    OFASTout = joinpath(res_path,"outD_HG$i")   # OpenFAST output file
    cont = parse_fxw(OFASTout, 5)

    #############################################################################################
    # Assign Variables
    NoRows = size(cont)[2]

    t, Surge, Heave, Pitch, Wave1Elev, Wave1Elev1, FAIRTEN1, FAIRTEN2, 
    HydroFxi, HydroFzi, HydroMyi, PtfmTAxt, PtfmTAzt = [cont[:,i] for i = 1:NoRows]

    aCoM = sqrt.(PtfmTAxt.^2 .+ PtfmTAzt.^2)  # CoM acceleration
    Nₜ = length(t)       
    dt = t[2]-t[1]

    #############################################################################################
    # Critical Event statistics & selection
    RespVar = FAIRTEN2[:]    # The response variable to be analysed
    case_str = "MaxFair"
    varname = "Fairlead tension"

    # # Response Statistics
    # StatMetrics = sample_stats(RespVar)    
    # push!(StatMetrics, rms(RespVar))
    # RefVar = StatMetrics[end]   # Reference/normalisation variable (rms, median, mean, RespVar[1])

    # Find maximum response value
    Fₘₐₓ[i] = findmax(RespVar)[1]
    iₘₐₓ = findmax(RespVar)[2]
    tₘₐₓ[i] = t[iₘₐₓ]

    # ## Wave Spectra
    # FR, MAG, _ = one_side_asp(Wave1Elev, t)               # OpenFAST - spectrum of surface elevation
    # FR_1st, MAG_1st, _ = one_side_asp(Wave1Elev1, t)      # OpenFAST - spectrum of 1st order elevation
    # FR_2nd, MAG_2nd, _ = one_side_asp(Wave1Elev .- Wave1Elev1, t)   # OpenFAST - spectrum 2nd order
    # fₚ = FR[findmax(MAG)[2]]

    plt_Wave1Elev = plot(xlab =L"t~[s]", ylab=L"\eta~[m]", title = "$(heads[5])", legend=:topleft, palette=[cb[11]])
    plot!(t, cont[:,5], lw=2, lab = "HyperGroup")

    plt_FAIRTEN2 = plot(xlab =L"t~[s]", ylab=L"F~[kN]", title = "$(heads[8])", legend=:topleft, palette=[cb[11]])
    plot!(t, cont[:,8], lw=2, lab = "HyperGroup")

    # ## Plot OpenFAST surface elevation spectra
    # plt_Spectra = plot(xlab = L"f~[Hz]", ylab = L"S(f)~[m]", palette=[cb[11];cb[4];cb[8]]) 
    # plot!(FR, MAG, lw=3, lab = "OpenFAST")
    # plot!(FR_1st, MAG_1st, lw=2, line =:dot, lab = L"1^{st}~\mathrm{order}")
    # plot!(FR_2nd, MAG_2nd, lw=2, line =:dot, lab = L"2^{nd}~\mathrm{order}")
    # plot!(xlim=(0,fcut),  minorgrid=true)

    plt_WF = plot(plt_Wave1Elev, plt_FAIRTEN2, layout = @layout[a; b])
    push!(AllPlots, plt_WF);    push!(plt_img_names, "Wave_Fair_$i")

    # ## Find all maxima above a given threshold of the selected response time history
    # RespVar_Peaks, _ = peaks_max_ext(RespVar.-RefVar,t, 0, 0, false)

    # ## Find and sort extreme maxima in the selected response time history
    # MinPeakVal = StatMetrics[1]    # Amplitude threshold
    # Aᵢ, tᵢ, i⁺, _ = peaks_max_ext(RespVar,t, MinPeakVal, 0, true)
    # Nce = length(Aᵢ)        # Number of peaks above threshold

    # ## Distribution of normalised peaks with reference to the rms value - Weibull
    # dataset = RespVar_Peaks./RefVar;     dist_type = Weibull;   bins = 2^9
    # pdf_values, fitted_dist, params, samp_params, hst_plt, QQplt, value_range = dist_fit(dataset, dist_type, bins)

    # kʷ, λʷ = params[1], params[2]
    # ptile = 1/100;  F̅₁ₚ = λʷ * (-log(ptile))^1/kʷ;     F₁ₚ = (F̅₁ₚ+1)*RefVar
end

Fₘₐₓ, iₒᵤₜ = remove_outliers(Fₘₐₓ,3)

box_F = boxplot(Fₘₐₓ, ylab=L"F_{max}")
push!(AllPlots, box_F);    push!(plt_img_names, "Fmax_box")

plt_scatter = plot(Fₘₐₓ, seriestype=:scatter, xlab="HG", ylab=L"F_{max}")
push!(AllPlots, plt_scatter);    push!(plt_img_names, "Fmax_scatter")

## Histogram of peak values
dataset = Fₘₐₓ;     dist_type = Normal;   bins = 2^9
pdf_values, fitted_dist, params, samp_params, hst_plt, QQplt, value_range = dist_fit(dataset, dist_type, bins)

plot!(hst_plt, title="PDF - P(Fₘₐₓ)", xlab=L"F_{max}", ylab=L"P(F_{max})")
FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])
push!(AllPlots, FitPlot);    push!(plt_img_names, "Fmax_hist")


# ### Plot PDF - Compare sample & fit
# hplt_F2 = histogram(SampVar, normalize=:pdf, lab=:false, palette=[cb[8]; cb[11]; cb[2]; cb[5]])
# plot!(hplt_F2, title="PDF - $varname maxima", xlab=L"\frac{F_p}{F_{rms}}-1", ylab="P(f)")
# plot!(hplt_F2, WeibRange, WeibPDF, lw=2, lab="Weibull: κ=$(round(kʷ*1e3)/1e3), λ=$(round(λʷ*1e3)/1e3)")
# plot!(hplt_F2,[mean(SampVar)+3*std(SampVar); mean(SampVar)+3*std(SampVar)], [0; maximum(WeibPDF)], ls=:dashdot, lab=L"\mu+3\sigma")
# plot!([F̅₁ₚ; F̅₁ₚ], [0; maximum(WeibPDF)], ls=:dashdot, lab="Q($(ptile*100)%)")
# plot!(twinx(),WeibRange, 1 .-WeibCDF, lw=2, ylab="Quantile", lab="Q(p)", ylim=(0,1), legendposition=:bottomright, palette=[cb[5]], ls=:dash)

#############################################################################################
# PLOTS
## Plot OpenFast results
for i ∈ 1:HG_fcount+1
    display(AllPlots[i])
end

display(AllPlots[end])
savefig(joinpath(res_path,"Fmax_hist.svg"))
savefig(joinpath(res_path,"Fmax_hist.png"))

end
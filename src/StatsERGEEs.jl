# Load necessary modules
using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles
using FFTW, DSP, Statistics, BSplineKit, StatsBase, Distributions, Copulas, Optim
using CurveFit: linear_fit
using StatsPlots
import ColorSchemes.darkrainbow

gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))
cb = darkrainbow

# Include necessary scripts for functions
include("func/text_process.jl"); include("func/signal_processing.jl")
include("func/peak_detect.jl");  include("func/runge_kutta.jl")
include("func/fugees.jl");       include("func/directories.jl")
include("func/wave_theory.jl");  include("func/jonswap.jl")
include("func/stat_process.jl")

#############################################################################################
# INPUT
const ρ, g, d::Float64 = 1025.0, 9.81, 100.0  # density [kg/m³], gravity [m/s²], depth [m]
const T₀ₛ, T₀ₚ, T₀ₕ::Float64 = 113.764, 26.253, 17.354 # FOWT eigenfrequencies (surge, pitch, heave)

# TO RUN ISOLATED CASE OF Hₛ,Tₚ, RESTART THE LANGUAGE SERVER (Alt+J > Alt+R), i.e. CLEAN MAIN SCOPE
Hₛⱼ, Tₚⱼ::Float64 = 5.0, 8.0    # Significant wave height [m], Peak wave period [m]
Tₑ, Tcut::Float64 = 32.0, 1.88  # Return and cut-off period [s] (lower and upper bounds of freq)
fcut = 1/Tcut   # Cut-off frequency

## Time vector (non-dimensional)
t̅ₑ = 2^9/Tₚⱼ        # Duration
Nₜ::Int64 = 2^12    # No of time steps

dt = t̅ₑ/(Nₜ-1)  # Time step
t̅ = zeros(Float64,Nₜ)   # Initialise non-dimensional time vector
[t̅[i] = (i-1)*dt-t̅ₑ/2 for i ∈ 1:Nₜ] # Non-dimensional time vector

# FLAGS
wout = true     # Write output files
wtot = false    # Write total stats - Use to store the summarised datasets from multiple cases

# Auxilliary variable to check whether statistics should be concatenated
prior_case_id = ""
if @isdefined(case_id)
    prior_case_id = case_id
end

# PATHS
case_id = js_case(Hₛⱼ, Tₚⱼ, 3.3, true)
libpath = joinpath(pwd(),"library")
stats_path = joinpath(libpath,"SE",case_id,"0")
case_str = "MaxFair"

# Decide if it is a multiple simulation (concatenate stats)
if prior_case_id !== case_id && prior_case_id !== ""
    global multisim = true
    stats_path = joinpath(libpath,"SE","TotStats")
else
    global multisim = false
end

#############################################################################################
# Read from OpenFAST simulation results
FirstRun, NoRuns = 1, 20
caseFairTen, Hₛᴿ, Tₚᴿ, ϵₛₛ, H̃ₛᴿ, T̃ₚᴿ, NoEv, TotNoEv, plt_CStats = read_ofast_stats(case_str,case_id,FirstRun,NoRuns)

if wout && !multisim
    col_head = ["median(Hₛ) [m]"; "median(Tₚ) [s]"; "median(ϵₛₛ)"; "median(Fₘₒₒᵣ) [kN]"; "TotNoEv"]
    samp_stats = [H̃ₛᴿ; T̃ₚᴿ; median(ϵₛₛ); median(caseFairTen); TotNoEv]
    fid = joinpath(stats_path,"OFAST_stats")
    open(fid, "w")
    writedlm(fid, [col_head samp_stats], '\t')
end

# Create box plots
if !@isdefined(box_Hs)
    global box_Hs, box_Tp, box_steep, box_FairTen
    box_Hs = plot(ylabel=L"H_s~[m]")
    box_Tp = plot(ylabel=L"T_p~[s]")
    box_steep = plot(ylabel=L"\epsilon_{ss}")
    box_FairTen = plot(ylabel=L"F_{moor}~[kN]")
end
boxplot!(box_Hs, Hₛᴿ, label=L"H^J_s="*"$(Hₛⱼ), "*L"T^J_p="*"$Tₚⱼ", legend=:outerright, legendcolumns=1)
boxplot!(box_Tp, Tₚᴿ, label=L"H^J_s="*"$(Hₛⱼ), "*L"T^J_p="*"$Tₚⱼ", legend=:outerright, legendcolumns=1)
boxplot!(box_steep, ϵₛₛ, label=L"H^J_s="*"$(Hₛⱼ), "*L"T^J_p="*"$Tₚⱼ", legend=:outerright, legendcolumns=1)
boxplot!(box_FairTen, caseFairTen, label=L"H^J_s="*"$(Hₛⱼ), "*L"T^J_p="*"$Tₚⱼ", legend=:outerright, legendcolumns=1)

# Dimensionless fairlead tension
# FT_thres = floor(minimum(caseFairTen))  
FT_thres = 611.389 # Initial fairlead tension [kN]
caseFairTen = caseFairTen./FT_thres

#############################################################################################
# Read and assign event and Gaussian Elementary Envelopes' (GEEs) parameters
# CWE_datasets, GEE_datasets, allΔtᶜ, allΔtᶜᵣ, Tᵖ₂₋,_,_,_,_,pltGEE5,_ = read_ergees_pars(case_str,case_id, FirstRun, NoRuns, TotNoEv, t̅)
CWE_case_dsets, GEE_case_dsets, Δt_case_dsets, other_datasets, Tᵖ₂₋,_,_,_,_,_,pltGEE5,_ = read_ergees_pars(case_str,case_id, FirstRun, NoRuns, TotNoEv, t̅)

CaseNoGEEs = size(GEE_case_dsets)[1] 
Øⁿ = Array{Float64}(undef,0,0)
# Initialise necessary variables if the main scope is clean
if !@isdefined(CWE_datasets)
    global CWE_datasets,GEE_datasets,Δt_datasets,FairTen
    CWE_datasets,GEE_datasets,Δt_datasets,FairTen = [Øⁿ[:,:] for _ = 1:4]
end
# Concatenate datasets in the case of combined simulations
if multisim
    CWE_datasets = vcat(CWE_datasets,hcat(CWE_case_dsets[:,:],caseFairTen[:]))
    GEE_datasets = vcat(GEE_datasets,GEE_case_dsets)
    Δt_datasets = vcat(Δt_datasets,Δt_case_dsets)
    FairTen = vcat(FairTen,caseFairTen)
else
    CWE_datasets = hcat(CWE_case_dsets[:,:],caseFairTen[:])
    GEE_datasets = GEE_case_dsets[:,:]
    Δt_datasets = Δt_case_dsets[:,:]
    FairTen = caseFairTen[:]
end

allAₒ, allTₒ, allTₒᵣ, allTₒₑ, alltᶜ, alltᶜᵣ, alltᶜₑ, allβ = [GEE_datasets[:,i] for i ∈ 1:8]
allΩ, allβ̇, allΔTₑᵥ, allN = [CWE_datasets[:,i] for i ∈ 1:4]
allΔtᶜ, allΔtᶜᵣ, allΔtᶜₑ = [Δt_datasets[:,i] for i ∈ 1:3]

# Create box plots
if !@isdefined(box_A)
    global box_A, box_T, box_Tᵣ, box_Ω, box_β̇, box_ΔTₑᵥ, box_N, box_Δtᶜ, box_Δtᶜᵣ
    box_A = plot(ylabel=L"\overline{A}")
    box_T = plot(ylabel=L"\overline{T}")
    box_Tᵣ = plot(ylabel=L"\overline{T}_d")
    box_Ω = plot(ylabel=L"\Omega")
    box_β̇ = plot(ylabel=L"\dot{\beta}")
    box_ΔTₑᵥ = plot(ylabel=L"\overline{\Delta T}_{ev}")
    box_N = plot(ylabel=L"\overline{N}")
    box_Δtᶜ = plot(ylabel=L"\overline{\Delta t}^c")
    box_Δtᶜᵣ = plot(ylabel=L"\overline{\Delta t}^c_d")
end
boxplot!(box_A, GEE_case_dsets[:,1], label=L"H^J_s="*"$(Hₛⱼ), "*L"T^J_p="*"$Tₚⱼ", legend=:outerright, legendcolumns=1)
boxplot!(box_T, GEE_case_dsets[:,2], label=L"H^J_s="*"$(Hₛⱼ), "*L"T^J_p="*"$Tₚⱼ", legend=:outerright, legendcolumns=1)
boxplot!(box_Tᵣ, GEE_case_dsets[:,3], label=L"H^J_s="*"$(Hₛⱼ), "*L"T^J_p="*"$Tₚⱼ", legend=:outerright, legendcolumns=1)
boxplot!(box_Ω, CWE_case_dsets[:,1], label=L"H^J_s="*"$(Hₛⱼ), "*L"T^J_p="*"$Tₚⱼ", legend=:outerright, legendcolumns=1)
boxplot!(box_β̇, CWE_case_dsets[:,2], label=L"H^J_s="*"$(Hₛⱼ), "*L"T^J_p="*"$Tₚⱼ", legend=:outerright, legendcolumns=1)
boxplot!(box_ΔTₑᵥ, CWE_case_dsets[:,3], label=L"H^J_s="*"$(Hₛⱼ), "*L"T^J_p="*"$Tₚⱼ", legend=:outerright, legendcolumns=1)
boxplot!(box_N, CWE_case_dsets[:,4], label=L"H^J_s="*"$(Hₛⱼ), "*L"T^J_p="*"$Tₚⱼ", legend=:outerright, legendcolumns=1)
boxplot!(box_Δtᶜ, Δt_case_dsets[:,1], label=L"H^J_s="*"$(Hₛⱼ), "*L"T^J_p="*"$Tₚⱼ", legend=:outerright, legendcolumns=1)
boxplot!(box_Δtᶜᵣ, Δt_case_dsets[:,2], label=L"H^J_s="*"$(Hₛⱼ), "*L"T^J_p="*"$Tₚⱼ", legend=:outerright, legendcolumns=1)
box_t = boxplot(GEE_case_dsets[:,5:7], label=[L"\overline{t^c}" L"\overline{t^c}_d" L"\overline{t^c}_{ev}"], ylab="Dimensionless time")

AllBoxPlots = [box_Hs box_Tp box_steep box_FairTen box_A box_T box_Tᵣ box_Ω box_β̇ box_ΔTₑᵥ box_N box_Δtᶜ box_Δtᶜᵣ box_t]
BoxPlotNames = ["box_Hs" "box_Tp" "box_steep" "box_Fair" "box_A" "box_T" "box_Tr" "box_Ω" "box_betadot" "box_DTev" "box_N" "box_Dtc" "box_Dtcr" "box_t"]

for i ∈ AllBoxPlots
    display(i)
end
savefig(box_t,joinpath(stats_path,"box_t.png"))
savefig(box_t,joinpath(stats_path,"box_t.svg"))
if wtot
    for i ∈ 1:length(AllBoxPlots)
        savefig(AllBoxPlots[i],joinpath(stats_path,"$(BoxPlotNames[i]).png"))
        savefig(AllBoxPlots[i],joinpath(stats_path,"$(BoxPlotNames[i]).svg"))
    end
end

#############################################################################################
# STATISTICAL ANALYSIS
## Statistics of Samples
### For critical wave event parameters (FairTen, Ω, ωₚ, β̇, T₂₋, ΔTₑᵥ, N̅, maxAₒ)
CWE_stats = sample_stats_multi(CWE_datasets)

### For GEEs' parameters (A̅ₒ,T̅ₒ,t̅ᶜₒ, β)
GEE_stats = sample_stats_multi(GEE_datasets)

### For Δt̅ᶜ datasets
Δt_stats = sample_stats_multi(Δt_datasets)

## Summarise sample statistics & write into text file
sample_statistics = round.(1000*[CWE_stats; GEE_stats; Δt_stats])./1000
chead = ["Variable" " x̅ " " s " "Med" "IQR" "Skew" "Kurt"]   # Column header
rhead = ["Ω"; "β̇";"ΔT̅ₑᵥ";"N̅";"max(A)";"FairTen";"A̅";"T̅";"T̅ᵣ";"T̅ₑ";"t̅ᶜ";"t̅ᶜᵣ";"t̅ᶜₑ";"β";"Δt̅ᶜ";"Δt̅ᶜᵣ";"Δt̅ᶜₑ"]  # Row header

if wout
    fid = joinpath(stats_path,"sample_stats")
    open(fid, "w")
    writedlm(fid, [chead; [rhead sample_statistics]], '\t')
end

# Store summarised datasets
if wtot && multisim && wout
    println("WARNING: Attempting to store total statistics.")
    println("WARNING: Previous results will be overwritten.")
    println("WARNING: Are you sure you want to proceed? (y/n)")
    uinp = readline()
    if uinp == "n"
        error("Script execution aborted!")
    else
        fid = joinpath(stats_path,"summary_CWE_datasets")
        open(fid, "w")
        writedlm(fid, [["Ω" "β̇" "ΔT̅ₑᵥ" "N̅" "max(A)" "FairTen"]; CWE_datasets], '\t')

        fid = joinpath(stats_path,"summary_GEE_datasets")
        open(fid, "w")
        writedlm(fid, [["A̅" "T̅" "T̅ᵣ" "T̅ₑ" "t̅ᶜ" "t̅ᶜᵣ" "t̅ᶜₑ" "β"]; GEE_datasets], '\t')

        fid = joinpath(stats_path,"summary_Δt_datasets")
        open(fid, "w")
        writedlm(fid, [["Δt̅ᶜ" "Δt̅ᶜᵣ" "Δt̅ᶜₑ"]; Δt_datasets], '\t')
    end
end

# FITTING STATISTICAL DISTRIBUTIONS
## Variable labels
GEE_fit_vars = ("A̅","T̅","T̅ᵣ","T̅ₑ","t̅ᶜ","t̅ᶜᵣ","t̅ᶜₑ","β")
CWE_fit_vars = ("Ω","β̇","ΔT̅ₑᵥ","N̅","max(A̅)","FairTen")
Δt_fit_vars = ("Δtᶜ","Δtᶜᵣ","Δt̅ᶜₑ")
## Variable labels in LaTeX notation
GEE_fvars_ltx = (L"\overline{A}",L"\overline{T}",L"\overline{T}_r",L"\overline{T}_e",L"\overline{t}^c",L"\overline{t}^c_r",L"\overline{t}^c_e",L"\beta")
CWE_fvars_ltx = (L"\Omega", L"\dot{\beta}", L"\Delta\overline{T}_{ev}", L"\overline{N}", L"max(\overline{A})",L"\overline{F}_{moor}")
Δt_fvars_ltx = (L"\Delta\overline{t}^c",L"\Delta\overline{t}^c_r",L"\Delta\overline{t}^c_e")
## Filenames of figures to save
GEE_ffnames = ("A_dist","T_dist","Tr_dist","Te_dist","tc_dist","tcr_dist","tce_dist","beta_dist")
CWE_ffnames = ("Om_dist","betadot_dist","DTev_dist","N_dist","maxA_dist","FairTen_dist")
Δt_ffnames = ("Dtc_dist", "Dtcr_dist","Dtcre_dist")

GEE_dist_types = (Weibull, LogNormal, LogNormal, LogNormal, Cauchy, Cauchy, Normal, Cauchy)
CWE_dist_types = (LogNormal, Normal, LogNormal, LogNormal, Weibull, LogNormal)
Δt_dist_types = (LogNormal, LogNormal, LogNormal)

IndexFitVar = 0
NoFitVars = length(GEE_fit_vars) + length(CWE_fit_vars) + length(Δt_fit_vars)
Nbins = 2^9
Dist_stats = zeros(Float64,NoFitVars,8) # [Π₁ Π₂ μ σ med() modes() skew() kurt()]
Øᵃⁿʸ = Array{Any}(undef,0)
AllFitPlots, GEE_fit_dist, CWE_fit_dist, Δt_fit_dist = [Øᵃⁿʸ[:] for _ = 1:4]
DistTypes, plt_img_names = Array{String}(undef,0), Array{String}(undef,0)

for i ∈ 1:length(GEE_fit_vars)
    global IndexFitVar += 1
    dataset = GEE_datasets[:,i];   dist_type = GEE_dist_types[i];
    pdf_values, fit_dist, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
    plot!(hst_plt, title="PDF - P($(GEE_fit_vars[i]))", xlab=GEE_fvars_ltx[i], ylab="P("*GEE_fvars_ltx[i]*")")
    FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])
    push!(AllFitPlots, FitPlot);    push!(plt_img_names, GEE_ffnames[i])
    push!(DistTypes,string(dist_type))
    Dist_stats[IndexFitVar,:] = [dist_params' mean(fit_dist) std(fit_dist) median(fit_dist) modes(fit_dist) skewness(fit_dist) kurtosis(fit_dist)]
    push!(GEE_fit_dist, fit_dist)
end
Aₒ_Fit, Tₒ_Fit, Tₒᵣ_Fit, tᶜ_Fit, tᶜᵣ_Fit, tᶜₑ_Fit = [GEE_fit_dist[i] for i ∈ 1:6]

for i ∈ 1:length(CWE_fit_vars)
    global IndexFitVar += 1
    dataset = CWE_datasets[:,i];   dist_type = CWE_dist_types[i];
    pdf_values, fit_dist, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
    plot!(hst_plt, title="PDF - P($(CWE_fit_vars[i]))", xlab=CWE_fvars_ltx[i], ylab="P("*CWE_fvars_ltx[i]*")")
    FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])
    push!(AllFitPlots, FitPlot);    push!(plt_img_names, CWE_ffnames[i])
    push!(DistTypes,string(dist_type))
    Dist_stats[IndexFitVar,:] = [dist_params' mean(fit_dist) std(fit_dist) median(fit_dist) modes(fit_dist) skewness(fit_dist) kurtosis(fit_dist)]
    push!(CWE_fit_dist, fit_dist)
end
ΩFit, β̇_Fit, N̅_Fit_Fit, ΔTₑᵥ_Fit = [CWE_fit_dist[i] for i ∈ 1:4]

for i ∈ 1:length(Δt_fit_vars)
    global IndexFitVar += 1
    dataset = Δt_datasets[:,i];   dist_type = Δt_dist_types[i];
    pdf_values, fit_dist, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
    plot!(hst_plt, title="PDF - P($(Δt_fit_vars[i]))", xlab=Δt_fvars_ltx[i], ylab="P("*Δt_fvars_ltx[i]*")")
    FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])
    push!(AllFitPlots, FitPlot);    push!(plt_img_names, Δt_ffnames[i])
    push!(DistTypes,string(dist_type))
    Dist_stats[IndexFitVar,:] = [dist_params' mean(fit_dist) std(fit_dist) median(fit_dist) modes(fit_dist) skewness(fit_dist) kurtosis(fit_dist)]
    push!(Δt_fit_dist, fit_dist)
end
Δtᶜ_Fit, Δtᶜᵣ_Fit, Δtᶜₑ_Fit = [Δt_fit_dist[i] for i ∈ 1:3]

# # Student-t Fit
# initial_params = [mean(alltᶜᵣ); std(alltᶜᵣ); 5.0]
# lower_bounds = [mean(alltᶜᵣ)-3*std(alltᶜᵣ); 0.0; 2.0]
# upper_bounds = [mean(alltᶜᵣ)+3*std(alltᶜᵣ); 20.0; 500.0]
# plt_hist, fitted_params, fitted_dist = student_fit(alltᶜᵣ, initial_params, lower_bounds, upper_bounds)
# display(plt_hist)

## Pearsons correlation coefficient
cvar = cov(hcat(GEE_datasets[:,4],GEE_datasets[:,1]))
sigA = std(GEE_datasets[:,4])
sigB = std(GEE_datasets[:,1])
rho = cvar ./ (sigA * sigB)

chead_D = ["Variable" "Distribution" "Π₁" "Π₂" " μ " "√σ²" "Med" "Mode" "Skew" "Kurt"]   # Column header
rhead_D = ["A̅";"T̅";"T̅ᵣ";"T̅ₑ";"t̅ᶜ";"t̅ᶜᵣ";"t̅ᶜₑ";"β";"Ω"; "β̇";"ΔT̅ₑᵥ";"N̅";"maxA";"FT";"Δt̅ᶜ";"Δt̅ᶜᵣ";"Δt̅ᶜₑ"]  # Row header

if wout
    fid = joinpath(stats_path,"Dist_stats")
    open(fid, "w")
    writedlm(fid, [chead_D; [rhead_D DistTypes round.(1000*Dist_stats)./1000]], '\t')
end

#############################################################################################
# Copulas
## C[F(Tₒ),F(Aₒ)]
JPD_TᵣA, Cop_TᵣA, plt_cop_TᵣA, plt_cop_TᵣA_OG = copula2d_fit_eval(GaussianCopula,allTₒᵣ,allAₒ,Tₒᵣ_Fit,Aₒ_Fit,2^16,100)
### Add labels to the copula density plots
plot!(plt_cop_TᵣA_OG, xlabel=L"\overline{T}_r", ylabel=L"\overline{A}")
plot!(plt_cop_TᵣA, xlabel=L"P(\overline{T}_r)", ylabel=L"P(\overline{A})")

## C[F(tᶜ),F(Aₒ)]
JPD_tᶜᵣA, Cop_tᶜᵣA, plt_cop_tᶜᵣA, plt_cop_tᶜᵣA_OG = copula2d_fit_eval(GaussianCopula,alltᶜᵣ,allAₒ,tᶜᵣ_Fit,Aₒ_Fit,2^16,100)
### Add labels to the copula density plots
plot!(plt_cop_tᶜᵣA_OG, xlabel=L"\overline{t}^c_r", ylabel=L"\overline{A}")
plot!(plt_cop_tᶜᵣA, xlabel=L"P(\overline{t}^c_r)", ylabel=L"P(\overline{A})")

## 3-variate copula: C[F(tᶜ), F(Tₒ), F(Aₒ)]
JPD_tᶜᵣTA, Cop_tᶜᵣTA, plt_3vc_tᶜTₒ, plt_3vc_tᶜTₒ_OG, plt_3vc_tᶜAₒ, plt_3vc_tᶜAₒ_OG, plt_3vc_TₒAₒ, plt_3vc_TₒAₒ_OG = copula3d_fit_eval(GaussianCopula,alltᶜᵣ,allTₒ, allAₒ,tᶜᵣ_Fit,Tₒ_Fit, Aₒ_Fit,2^16,100)
### Add labels to the copula density plots
plot!(plt_3vc_tᶜTₒ_OG, xlabel=L"\overline{t}^c", ylabel=L"\overline{T}_o")
plot!(plt_3vc_tᶜAₒ_OG, xlabel=L"\overline{t}^c", ylabel=L"\overline{A}_o")
plot!(plt_3vc_TₒAₒ_OG, xlabel=L"\overline{T}_o", ylabel=L"\overline{A}_o")
plot!(plt_3vc_tᶜTₒ, xlabel=L"P(\overline{t}^c)", ylabel=L"P(\overline{T}_o)")
plot!(plt_3vc_tᶜAₒ, xlabel=L"P(\overline{t}^c)", ylabel=L"P(\overline{A}_o)")
plot!(plt_3vc_TₒAₒ, xlabel=L"P(\overline{T}_o)", ylabel=L"P(\overline{A}_o)")

plt_cop_tᶜᵣTA_OG = plot(plt_3vc_tᶜTₒ_OG, plt_3vc_tᶜAₒ_OG, plt_3vc_TₒAₒ_OG, layout= @layout[a b; c d], title = "3-variate copula f(tᶜ,T,A)")
plt_cop_tᶜᵣTA = plot(plt_3vc_tᶜTₒ, plt_3vc_tᶜAₒ, plt_3vc_TₒAₒ, layout= @layout[a b; c d], title = "3-variate copula f(tᶜ,T,A)")

#############################################################################################
### Sample from resulting joint probability distribution P(x₁,x₂)
# Binning of of the conditioning variable x₁
## Range for x₁ (lower and upper bound)   
x₁_hrange = mean(alltᶜᵣ)+std(alltᶜᵣ)   

Δtᶜ_μ = Dist_stats[end-1,1]
Δtᶜ_σ = Dist_stats[end-1,2]
Δtᶜ_var = (exp(Δtᶜ_σ^2)-1)*exp(2*Δtᶜ_μ+Δtᶜ_σ^2)
σᴮ = sqrt(Δtᶜ_var)    # Use of standard deviation of appropriate dataset to determine width of bin
wᴮ = 2*σᴮ   # Width of bin

x₁ˡᵇ = -(Int(floor(x₁_hrange/wᴮ))+0.5)*wᴮ   
x₁ᵘᵇ = (Int(floor(x₁_hrange/wᴮ))+0.5)*wᴮ 
Nᴮ = Int(floor((x₁ᵘᵇ-x₁ˡᵇ)/wᴮ)) + 1  # Number of bins
x₁_bins = range(x₁ˡᵇ,x₁ᵘᵇ,Nᴮ)   # Central values of bins for the x₁ range

plt_tcbins = plot(xlab=L"\overline{t}^c_r", ylab=L"\overline{A}")
for i ∈ 1:Nᴮ
    plot!(plt_tcbins, [x₁_bins[i]; x₁_bins[i]], [0; maximum(allAₒ)], ls=:dash, leg=:false)
end

# Compute the distribution in each tᶜ bin based on the A samples that belong to it
Aₒ_WeiFits = zeros(Float64,Nᴮ,2)
Aₒ_WeiModes = zeros(Float64,Nᴮ)
# plt_BinDist = plot(xlab=L"t^c", ylab=L"\overline{A}_o", zlab=L"P(\overline{A}_o)")
plt_BinDist = plot(xlab=L"\overline{A}", ylab=L"P(\overline{A})")
# Distribution of Aₒ per tᶜ bin
for j ∈ 1:Nᴮ
    x₁ᶜ = x₁ˡᵇ + j*wᴮ - wᴮ/2 # Which bin?
    tᶜ_samp, Aₒ_samp = [], []
    for i ∈ 1:length(allAₒ)
        if alltᶜ[i] < (x₁ᶜ+wᴮ/2) && alltᶜ[i] > (x₁ᶜ-wᴮ/2)
            push!(tᶜ_samp,alltᶜ[i])
            push!(Aₒ_samp,allAₒ[i])
        end
    end
    tᶜ_samp, Aₒ_samp = Float64.(tᶜ_samp), Float64.(Aₒ_samp)
    println("Number of samples in this bin: $(length(Aₒ_samp))")

    dataset = Aₒ_samp;   dist_type = Weibull;  Nbins = 2^9
    pdf_values, fitted_dist, dist_params, samp_params, hst_plt, QQplt, Aₒ_Wrang = dist_fit(dataset, dist_type, Nbins)
    # plot3d!(plt_BinDist, x₁ᶜ*ones(Nbins), Aₒ_Wrang, pdf_values, lw=2, lab="t= $(round(x₁ᶜ*1000)/1000)")
    plot!(plt_BinDist, Aₒ_Wrang, pdf_values, lw=2, lab="t= $(round(x₁ᶜ*1000)/1000)")
    # plot3d!(hst_plt,  title="PDF - P(A̅ₒ|tᶜ bin = $(round(x₁ᶜ*1000)/1000)", xlab=L"\overline{A}_o", ylab=L"P(\overline{A}_o)")
    # display(hst_plt)

    Aₒ_WeiFits[j,:] = dist_params[:]
    kʷ = dist_params[1]
    λʷ = dist_params[2]
    mode_value = λʷ*((kʷ-1)/kʷ)^(1/kʷ)
    Aₒ_WeiModes[j] = mode_value
    # println("Weibull Mode Value: $mode_value")
end

# Store final A Weibull fits per tᶜ bin
if wtot && multisim && wout
    fid = joinpath(stats_path,"A_Weibull_Fits")
    open(fid, "w")
    writedlm(fid, [["t̅ᵇ" "κ" "λ" "mode"]; [x₁_bins Aₒ_WeiFits Aₒ_WeiModes]], '\t')
end

# Conditional sampling from the tᶜ-A joint pdf for a value of tᶜ
x₁ᶜ = x₁_bins[Int(floor(Nᴮ/2))+1] # Which bin?
x₁_samp, x₂_samp, x₁_range, hist_x2, hist_x1x2 = conditional_sampling(JPD_tᶜᵣA,10000,x₁ᶜ,wᴮ)

dataset = x₂_samp;   dist_type = Weibull;  Nbins = 2^9
pdf_values, fitted_dist, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
plot!(hst_plt,  title="PDF - P(A̅ₒ|tᶜ bin = $(round(x₁ᶜ*1000)/1000)", xlab=L"\overline{A}_o", ylab=L"P(\overline{A}_o)")
FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])

#############################################################################################
# Histograms
tcrbins = range(-t̅ₑ*Tₚⱼ/T₀ₚ, t̅ₑ*Tₚⱼ/T₀ₚ, length=154)
tcebins = range(-1, 1, length=25)
Aobins = range(0, 1.2, length=25)
Tobins = range(0, 3.5, length=35)

### Tᵣ-A
hplt_TA = histogram2d(allTₒᵣ, allAₒ, normalize=:pdf, bins=:sqrt, show_empty_bins=true, color=:plasma, xlab = L"\overline{T}_r", ylab = L"\overline{A}")

### tᶜᵣ-Aₒ
hist = fit(Histogram, (alltᶜᵣ, allAₒ), (tcrbins,Aobins), closed=:right)
etc = hist.edges[1]
eA = hist.edges[2]
probabs = (hist.weights / sum(hist.weights))'
hplt_tcA = heatmap(etc[1:end-1],eA[1:end-1],probabs, color=:plasma, xlab = L"\overline{t}_r", ylab = L"\overline{A}")

### tᶜₑ-Aₒ
hplt_tceA = histogram2d(alltᶜₑ, allAₒ, normalize=:pdf, bins=(tcebins,Aobins), show_empty_bins=true, color=:plasma, xlab = L"\overline{t}_e", ylab = L"\overline{A}")

### tᶜ-Tₒ
hplt_tcT = histogram2d(alltᶜᵣ, allTₒ, normalize=:pdf, bins=(tcrbins,Tobins), show_empty_bins=true, color=:plasma, xlab = L"\overline{t}_r", ylab = L"\overline{T}")

### tᶜₑ -Tₒ
hplt_tceT = histogram2d(alltᶜₑ, allTₒ, normalize=:pdf, bins=(tcebins,Tobins), show_empty_bins=true, color=:plasma, xlab = L"\overline{t}_e", ylab = L"\overline{T}")

#############################################################################################
#### Interpolation of most probable envelope shape from t-A heatmap
Aint = zeros(Float64, size(probabs)[2])

for i ∈ 1:size(probabs)[2]
    j = findmax(probabs[:,i])[2]
    Aint[i] = eA[j]
end

A_itp = interpolate(etc[1:end-1], Aint, BSplineOrder(4))    # B-spline interpolation
Aᵢ = A_itp.(t̅)  # Interpolation for entire t̅ range
# Truncate interpolation for subset of t̅ range
Aₒᵢ = zeros(Float64,Nₜ)
Aₒᵢ[Int(Nₜ/2-2^10):Int(Nₜ/2+2^10)] = Aᵢ[Int(Nₜ/2-2^10):Int(Nₜ/2+2^10)]

#############################################################################################
# Linear fit of all β angles
fitβ̃ = linear_fit(alltᶜ,allβ)
fitβ = linear_fit(alltᶜ.*Tₚⱼ,allβ)

# β̃reg = fitβ̃[1] .+ fitβ̃[2].*alltᶜ
β̃reg = fitβ̃[2].*alltᶜ
plot!(pltGEE5,alltᶜ,β̃reg)

#############################################################################################
# PLOTS
for i ∈ 1:length(AllFitPlots)
    display(AllFitPlots[i])
    savefig(AllFitPlots[i],joinpath(stats_path,"$(plt_img_names[i]).png"))
    savefig(AllFitPlots[i],joinpath(stats_path,"$(plt_img_names[i]).svg"))
end

display(plt_tcbins)
savefig(plt_tcbins, joinpath(stats_path,"tc_bins.png"))
savefig(plt_tcbins, joinpath(stats_path,"tc_bins.svg"))
display(plt_BinDist)
savefig(plt_BinDist, joinpath(stats_path,"A_bin_dists.png"))
savefig(plt_BinDist, joinpath(stats_path,"A_bin_dists.svg"))

# display(pltΩ)
# savefig(pltΩ, joinpath(stats_path,"om_tilde.png"))
# savefig(pltΩ, joinpath(stats_path,"om_tilde.svg"))

# display(pltT2)
# savefig(pltT2, joinpath(stats_path,"T2-.png"))
# savefig(pltT2, joinpath(stats_path,"T2-.svg"))

# display(pltΔT)
# savefig(pltΔT, joinpath(stats_path,"DTev.png"))
# savefig(pltΔT, joinpath(stats_path,"DTev.svg"))

# display(pltFair)
# savefig(pltFair, joinpath(stats_path,"FTstats.png"))
# savefig(pltFair, joinpath(stats_path,"FTstats.svg"))

display(plt_CStats)
savefig(plt_CStats, joinpath(stats_path,"fairten_stats.png"))
savefig(plt_CStats, joinpath(stats_path,"fairten_stats.svg"))
# display(pltGEE1)
# savefig(pltGEE1, joinpath(stats_path,"EWGpars3D.png"))   
# savefig(pltGEE1, joinpath(stats_path,"EWGpars3D.svg"))
# display(pltGEE2)
# savefig(pltGEE2, joinpath(stats_path,"To-Ao.png"))
# savefig(pltGEE2, joinpath(stats_path,"To-Ao.svg"))

# display(pltGEE3)
# savefig(pltGEE3, joinpath(stats_path,"tc-To.png"))
# savefig(pltGEE3, joinpath(stats_path,"tc-To.svg"))
# display(pltGEE4)
# savefig(pltGEE4, joinpath(stats_path,"tc-Ao.png"))
# savefig(pltGEE4, joinpath(stats_path,"tc-Ao.svg"))
display(pltGEE5)
savefig(pltGEE5, joinpath(stats_path,"t-beta.png"))
savefig(pltGEE5, joinpath(stats_path,"t-beta.svg"))
# display(pltGEE6)
# savefig(pltGEE6, joinpath(stats_path,"tc_diff.png"))
# savefig(pltGEE6, joinpath(stats_path,"tc_diff.svg"))

display(hplt_TA)
savefig(hplt_TA, joinpath(stats_path,"hist_To-Ao.png"))
savefig(hplt_TA, joinpath(stats_path,"hist_To-Ao.svg"))
display(hplt_tcA)
savefig(hplt_tcA, joinpath(stats_path,"hist_tc-Ao.png"))
savefig(hplt_tcA, joinpath(stats_path,"hist_tc-Ao.svg"))
display(hplt_tceA)
savefig(hplt_tceA, joinpath(stats_path,"hist_tce-Ao.png"))
savefig(hplt_tceA, joinpath(stats_path,"hist_tce-Ao.svg"))
display(hplt_tcT)
savefig(hplt_tcT, joinpath(stats_path,"hist_tc-To.png"))
savefig(hplt_tcT, joinpath(stats_path,"hist_tc-To.svg"))
display(hplt_tceT)
savefig(hplt_tceT, joinpath(stats_path,"hist_tce-To.png"))
savefig(hplt_tceT, joinpath(stats_path,"hist_tce-To.svg"))
display(plt_cop_TᵣA_OG)
savefig(plt_cop_TᵣA_OG, joinpath(stats_path,"Ao-To_copula.png"))
savefig(plt_cop_TᵣA_OG, joinpath(stats_path,"Ao-To_copula.svg"))
display(plt_cop_tᶜᵣA_OG)
savefig(plt_cop_tᶜᵣA_OG, joinpath(stats_path,"tc-Ao_copula.png"))
savefig(plt_cop_tᶜᵣA_OG, joinpath(stats_path,"tc-Ao_copula.svg"))
display(plt_cop_tᶜᵣTA_OG)
savefig(plt_cop_tᶜᵣTA_OG, joinpath(stats_path,"var3_copula.png"))
savefig(plt_cop_tᶜᵣTA_OG, joinpath(stats_path,"var3_copula.svg"))

# display(pltSP2)
# savefig(pltSP2, joinpath(stats_path,"sp_2nd-.png"))
# savefig(pltSP2, joinpath(stats_path,"sp_2nd-.svg"))
# display(pltPHI2)
# savefig(pltPHI2, joinpath(stats_path,"phi_2nd-.png"))
# savefig(pltPHI2, joinpath(stats_path,"phi_2nd-.svg"))
# display(plt_eta2)
# savefig(plt_eta2, joinpath(stats_path,"eta_2nd-.png"))
# savefig(plt_eta2, joinpath(stats_path,"eta_2nd-.svg"))

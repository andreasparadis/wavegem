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
const ρ, g, d::Float64 = 1025.0, 9.81, 100.0  # density [kg/m³], gravity [m/s²], depth [m]
Tₚⱼ = 8.0
Hₛⱼ = 5.0
Tₑ = 32
Tcut = 1.88
fcut = 1/Tcut

regress = false
wout = true

libpath = joinpath(pwd(),"library")
case_str = "MaxFair"
case_id = js_case(Hₛⱼ, Tₚⱼ, 3.3, true)

#############################################################################################
pltG = plot(xlab = L"\overline{t}", ylab = L"\overline{A}", palette=:darkrainbow)
pltGEE1 = plot3d(xlab = L"\overline{T}", ylab = L"\overline{t}_c", zlab = L"\overline{A}", legend=:false)
pltGEE2 = plot(xlab = L"\overline{T}", ylab = L"\overline{A}", legend=:false)
pltGEE3 = plot(xlab = L"\overline{t}_c", ylab = L"\overline{T}", legend=:false)
pltGEE4 = plot(xlab = L"\overline{t}_c", ylab = L"\overline{A}", legend=:false)
pltGEE5 = plot(xlab = L"\overline{t}_c", ylab = L"\tilde{\beta}"*" [rad]", legend=:false)
pltGEE6 = plot(xlab = "Event", ylab = L"\Delta\overline{t}_c", legend=:false)
pltSP2 = plot(xlab = "f [Hz]", ylab = L"S [m^2 s]", palette=:darkrainbow)
pltPHI2 = plot(xlab = "f [Hz]", ylab = L"\phi [rad]", palette=:darkrainbow)

t̅ₑ = 2^9/Tₚⱼ
Nₜ = 2^12
Nₛ₂ = Int(nextpow(2,Nₜ)/2)
dt = t̅ₑ/(Nₜ-1)
t̅ = zeros(Float64,Nₜ)
[t̅[i] = (i-1)*dt-t̅ₑ/2 for i ∈ 1:Nₜ]

FirstRun = 1
NoRuns = 20
NoEv = zeros(Int8,NoRuns)
TotNoEv, RunCnt = 0, 0
FairTen = []
CaseStats = zeros(Float64,NoRuns,6)
Hₛᴿ = zeros(Float64,NoRuns)
Tₚᴿ = zeros(Float64,NoRuns)

for run_id ∈ FirstRun:(FirstRun+NoRuns-1)
    global RunCnt += 1
    postOFpath = joinpath(libpath,"SE","$(case_id)","$(run_id)","0","postOFAST")
    fid_tinsts = joinpath(postOFpath,"MaxFair_tinsts_1ptile")
    # fid_tinsts = joinpath(postOFpath,"MaxFair_tinsts")
    Decpath = joinpath(libpath,"SE","$(case_id)","$(run_id)","Decomposition")
    fsimpar = joinpath(Decpath,"sim_pars")       # Sea-state parameters file path
    instances = parse_fxw(fid_tinsts, 0)
    fid_mfstat = joinpath(postOFpath,case_str,"StatMetrics")
    CaseStats[RunCnt,:] = parse_fxw(fid_mfstat, 1) # CaseStats = [μ σ med() rms() skew() kurt()]
    SimPars = parse_fxw(fsimpar,1)
    Hₛᴿ[RunCnt] = SimPars[1]
    Tₚᴿ[RunCnt] = SimPars[2]

    NoEv[RunCnt] = size(instances)[1]
    global TotNoEv += NoEv[RunCnt]
    for i ∈ 1:NoEv[RunCnt]
        push!(FairTen, instances[i,2])
    end
end

H̃ₛᴿ = mean(Hₛᴿ)
T̃ₚᴿ = mean(Tₚᴿ)

# Candlestick plot for Fairled tensions per run
plt_CStats = plot(legend=false, title="Fairlead tension response", xlabel="Simulation", ylabel=L"\mu \pm \sigma [N]")
for i in 1:NoRuns
    plot!(plt_CStats, [i, i], [CaseStats[i,1], CaseStats[i,1]], color=:blue, linewidth=3) # mean
    plot!(plt_CStats, [i, i], [CaseStats[i,1] - CaseStats[i,2], CaseStats[i,1] + CaseStats[i,2]], color=:black, linewidth=1) # ±std
end
plot!(plt_CStats, ylim=(0, ylims(plt_CStats)[2]))

# Hₛ of runs 
plot(Hₛᴿ, seriestype=:scatter, lab=:false, title="Hs of simulations")
# Tₚ of runs
plot(Tₚᴿ, seriestype=:scatter, lab=:false, title="Tp of simulations")

#############################################################################################
MaxNoEv = maximum(NoEv)
allG = zeros(Float64,Nₜ,TotNoEv)
all_sp2, all_phi = zeros(Float64,Nₛ₂,TotNoEv), zeros(Float64,Nₛ₂,TotNoEv)
fsp2 = range(0.0,4.0,Nₛ₂)
allΩ, allωₚ, allβ, allβ̇, alltᶜ, allAₒ, allTₒ, allT₂₋, allΔTₑᵥ, allΔtᶜ, allN, maxAₒ  = [], [], [], [], [], [], [], [], [], [], [], []
TotEvID = 0

for run_id ∈ FirstRun:(FirstRun+NoRuns-1)
    Decpath = joinpath(libpath,"SE","$(case_id)","$(run_id)","Decomposition")
    for evID ∈ 1:NoEv[run_id-FirstRun+1]
        evdir = joinpath(case_str,"EV$evID")
        Eventpath = joinpath(Decpath,evdir)
        GRrespath = joinpath(Eventpath,"ERGEEs")
        fid_pEWG = joinpath(GRrespath,"EWG_norm_pars")
        fid_pG = joinpath(GRrespath,"G_pars")
        fid_evpars = joinpath(Eventpath,"ev_pars")
        fid_eta2nd = joinpath(Eventpath,"eta_2nd-")
        # fid_eta2nd = joinpath(Eventpath,"event_lin")

        # Read parameters of Gaussian Regression of wave event
        cont0 = parse_fxw(fid_evpars, 1)     # Read event parameters
        T₂₋, ΔTₑᵥ = cont0[4], cont0[5]
        # if T₂₋ ≥ 0 && T₂₋ ≤ 1024
        push!(allT₂₋,T₂₋)
        push!(allΔTₑᵥ, ΔTₑᵥ)
        global TotEvID += 1

        cont = parse_fxw(fid_pEWG, 1)   # Read EWGs parameters
        t̅ᶜ, T̅ₒ, A̅ₒ, β̃ = cont[:,1], cont[:,2], cont[:,3], cont[:,4]
        lenT = length(t̅ᶜ)
        cont2 = parse_fxw(fid_pG, 1)     # Read event and GR parameters
        Hₛ, Tₚ, Ω = cont2[1:3]
        linfit_β̃ = linear_fit(t̅ᶜ,β̃)

        if !isnan(linfit_β̃[2])
            push!(allβ̇,linfit_β̃[2])
        end
        push!(allΩ, Ω)
        push!(allωₚ, 2π/Tₚ)
        push!(maxAₒ, maximum(A̅ₒ))
        push!(allN, lenT/ΔTₑᵥ)
        

        for EWG ∈ 1:lenT
            push!(allAₒ,A̅ₒ[EWG])
            push!(allTₒ,T̅ₒ[EWG])
            push!(alltᶜ,t̅ᶜ[EWG])
            push!(allβ,β̃[EWG])
        end

        # Peak instance intervals
        s_t̅ᶜ = sort(t̅ᶜ)
        for EWG ∈ 1:lenT-1
            push!(allΔtᶜ,s_t̅ᶜ[EWG+1]-s_t̅ᶜ[EWG])
        end

        # Resulting Gaussian approximation of the envelope
        G̅ = gauss_fun(t̅, A̅ₒ, t̅ᶜ, T̅ₒ)
        allG[:,TotEvID] = G̅[:]

        # Read 2nd minus spectrum
        cont_2nd = parse_fxw(fid_eta2nd, 0)     
        time_ev = cont_2nd[:,1]
        eta2nd = cont_2nd[:,2]
        freq_2nd, mag_2nd,phi_2nd,_ = one_side_asp(eta2nd, time_ev)

        # Re-sampling spectrum
        itp_sp2 = interpolate(freq_2nd, mag_2nd, BSplineOrder(4))
        itp_mag_2nd = itp_sp2.(fsp2)
        # Re-sampling phase
        itp_phi = interpolate(freq_2nd, unwrap(phi_2nd), BSplineOrder(4))
        itp_phi_2nd = itp_phi.(fsp2)

        all_sp2[:, TotEvID] = itp_mag_2nd[:]
        all_phi[:, TotEvID] = itp_phi_2nd[:]

        plot!(pltG, t̅, G̅, lab=:false, line=:dot)
        plot!(pltGEE1, [T̅ₒ], [t̅ᶜ], [A̅ₒ], seriestype=:scatter)
        plot!(pltGEE2, [T̅ₒ], [A̅ₒ], seriestype=:scatter, xlim=(0,xlims(pltGEE1)[2]), ylim=(0,zlims(pltGEE1)[2]))
        plot!(pltGEE3, [t̅ᶜ], [T̅ₒ], seriestype=:scatter)
        plot!(pltGEE4, [t̅ᶜ], [A̅ₒ], seriestype=:scatter)
        plot!(pltGEE5, [t̅ᶜ], [β̃], seriestype=:scatter)
        # plot!(pltGEE5, [t̅ᶜ], [unwrap(β̃)], seriestype=:scatter)
        plot!(pltGEE6, TotEvID*ones(Int64,lenT-1), diff(sort(t̅ᶜ)), seriestype=:scatter)
        # plot!(pltGEE6, [TotEvID], [mean(diff(sort(t̅ᶜ)))], seriestype=:scatter)
        # plot!(pltGEE6, [ΔTₑᵥ], [std(diff(sort(t̅ᶜ)))], seriestype=:scatter)
        plot!(pltSP2, fsp2, itp_mag_2nd, lab=:false, line=:dot)
        plot!(pltPHI2, fsp2, itp_phi_2nd, lab=:false, line=:dot)
        # end
    end
end

Gmean = mean!(ones(Nₜ,1),allG[:,1:TotEvID])
plot!(pltG,t̅, Gmean, lw=3, lab="Mean")

sp2mean = mean!(ones(Nₛ₂,1),all_sp2[:,1:TotEvID])
sp2mean = sp2mean[:,1]
plot!(pltSP2, fsp2, sp2mean, lw=3, lab="Mean", xlim=(0,0.5))

phi2mean = mean!(ones(Nₛ₂,1),all_phi[:,1:TotEvID])
plot!(pltPHI2, fsp2, phi2mean, lw=3, lab="Mean")

Tᵖ₂₋ = 1/fsp2[findmax(sp2mean)[2]]

H2 = Nₛ₂*sp2mean.*exp.(1im*phi2mean)
# H2 = Nₛ₂*sp2mean.*exp.(1im*2π*rand(Nₛ₂))
H2 = vcat(H2,0,conj(H2[end:-1:2]))
η₂₋ = real(ifft(H2))
η₂₋ = [η₂₋[Int(Nₜ/2):end]; η₂₋[1:Int(Nₜ/2-1)]]

t = t̅*Tₚⱼ
# η = real(Hₛⱼ*Gmean.*exp.(-1im * mean(allΩ)*2π/Tₚⱼ * t̅*Tₚⱼ))

plt_eta2 = plot(t, η₂₋, xlab = "t [s]", ylab = L"\eta_{2nd} ~[m]", lw=2)

#############################################################################################
# STATISTICAL ANALYSIS
allΩ = Float64.(allΩ);  allβ = Float64.(allβ);  allβ̇ = Float64.(allβ̇);  allN = Float64.(allN)
allAₒ = Float64.(allAₒ);    allTₒ = Float64.(allTₒ);    alltᶜ = Float64.(alltᶜ)
allΔtᶜ = Float64.(allΔtᶜ);  allT₂₋ = Float64.(allT₂₋);  allΔTₑᵥ = Float64.(allΔTₑᵥ) 
FairTen = Float64.(FairTen)
GEE_stats = mean([allAₒ allTₒ alltᶜ], dims=1)
std([allAₒ allTₒ alltᶜ], dims=1)

## Scatter plots
pltΩ = plot([allΩ], seriestype=:scatter, legend=:false ,xlab="Event", ylab=L"\Omega=\frac{\tilde{\omega}}{\omega_p}~[-]")
pltT2 = plot([allT₂₋], seriestype=:scatter, legend=:false ,xlab="Event", ylab=L"T_{2^-}~[s]")
pltΔT = plot([allΔTₑᵥ], seriestype=:scatter, legend=:false ,xlab="Event", ylab=L"\Delta T_{ev}~[s]")
pltFair = plot([FairTen], seriestype=:scatter, legend=:false ,xlab="Event", ylab="Fairlead Tension [kN]")
pltN = plot([allN], seriestype=:scatter, legend=:false ,xlab="Event", ylab=L"\frac{No WGs}{\Delta T_{ev}}")

### FairTen
dataset = FairTen;   dist_type = LogNormal;  Nbins = 2^9
pdf_values, fitted_dist, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
plot!(hst_plt, title="PDF - P(FairTen)", xlab=L"F~[N]", ylab=L"P(F)")
FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])
AllFitPlots = [FitPlot]
plt_img_names = ["Fairten_dist"]

### Event Durations
dataset = allΔTₑᵥ;   dist_type = LogNormal;  Nbins = 2^9
pdf_values, fitted_dist, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
plot!(hst_plt, title="PDF - P(ΔTₑᵥ)", xlab=L"\Delta T_{ev}~[s]", ylab=L"P(\Delta T_{ev})")
FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])
push!(AllFitPlots, FitPlot);    push!(plt_img_names, "DTev_dist")

## Histograms
### Ω
dataset = allΩ;   dist_type = LogNormal;  Nbins = 2^9
pdf_values, ΩFit, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
plot!(hst_plt, title="PDF - P(Ω)", xlab=L"Ω", ylab=L"P(Ω)")
FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])
push!(AllFitPlots, FitPlot);    push!(plt_img_names, "om_tilde_dist")

### β̇
dataset = allβ̇ ;   dist_type = Normal;  Nbins = 2^9
pdf_values, fitted_dist, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
plot!(hst_plt, title="PDF - P(β̇)", xlab=L"\dot{\beta}~[rad]", ylab=L"P(\dot{\beta})")
FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])
push!(AllFitPlots, FitPlot);    push!(plt_img_names, "beta_dot_dist")

### N
dataset = allN;   dist_type = LogNormal;  Nbins = 2^9
pdf_values, fitted_dist, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
plot!(hst_plt, title="PDF: P(N̅)", xlab=L"\overline{N}", ylab=L"P(\overline{N})")
FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])
push!(AllFitPlots, FitPlot);    push!(plt_img_names, "NoWGs")

### Δtᶜ
dataset = allΔtᶜ;   dist_type = LogNormal;  Nbins = 2^9
pdf_values, Δtᶜ_Fit, Δtᶜ_dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
plot!(hst_plt, title="PDF - P(Δtᶜ)", xlab=L"\overline{\Delta t^c}", ylab=L"P(\overline{\Delta t^c})")
FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])
push!(AllFitPlots, FitPlot);    push!(plt_img_names, "tc_diff_dist")

### Tₒ
dataset = allTₒ;   dist_type = LogNormal;  Nbins = 2^9
pdf_values, Tₒ_Fit, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
plot!(hst_plt, title="PDF - P(T̅ₒ)", xlab=L"\overline{T}_o", ylab=L"P(\overline{T}_o)")
FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])
push!(AllFitPlots, FitPlot);    push!(plt_img_names, "To_dist")

### tᶜ
dataset = alltᶜ;   dist_type = Cauchy;  Nbins = 2^9
pdf_values, tᶜ_Fit, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
plot!(hst_plt,  title="PDF - P(t̅ᶜ)", xlab=L"\overline{t}^c", ylab=L"P(\overline{t}^c)")
FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])
push!(AllFitPlots, FitPlot);    push!(plt_img_names, "tc_dist")
# # Student-t Fit
# initial_params = [mean(alltᶜ); std(alltᶜ)/10; 5.0]
# lower_bounds = [mean(alltᶜ)-3*std(alltᶜ); 0.0; 2.0]
# upper_bounds = [mean(alltᶜ)+3*std(alltᶜ); 10.0; 500.0]
# plt_hist, fitted_params, fitted_dist = student_fit(alltᶜ, initial_params, lower_bounds, upper_bounds)
# display(plt_hist)

### Aₒ
dataset = allAₒ;   dist_type = Weibull;  Nbins = 2^9
pdf_values, Aₒ_Fit, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
plot!(hst_plt,  title="PDF - P(A̅ₒ)", xlab=L"\overline{A}_o", ylab=L"P(\overline{A}_o)")
FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])
push!(AllFitPlots, FitPlot);    push!(plt_img_names, "Ao_dist")

cvar = cov(hcat(alltᶜ,allAₒ))
sigA = std(alltᶜ)
sigB = std(allAₒ)

## Pearson correlation coefficient
rho = cvar ./ (sigA * sigB)
#############################################################################################
# Copulas
## C[F(Tₒ),F(Aₒ)]
JPD_TₒAₒ, Cop_TₒAₒ, plt_cop_TₒAₒ, plt_cop_TₒAₒ_OG = copula2d_fit_eval(GaussianCopula,allTₒ,allAₒ,Tₒ_Fit,Aₒ_Fit,2^16,100)
### Add labels to the copula density plots
plot!(plt_cop_TₒAₒ_OG, xlabel=L"\overline{T}_o", ylabel=L"\overline{A}_o")
plot!(plt_cop_TₒAₒ, xlabel=L"P(\overline{T}_o)", ylabel=L"P(\overline{A}_o)")

## C[F(Aₒ),F(Tₒ)]
JPD_AₒTₒ, Cop_AₒTₒ, plt_cop_AₒTₒ, plt_cop_AₒTₒ_OG = copula2d_fit_eval(GaussianCopula,allAₒ,allTₒ,Aₒ_Fit,Tₒ_Fit,2^16,100)
### Add labels to the copula density plots
plot!(plt_cop_AₒTₒ_OG, xlabel=L"\overline{A}_o", ylabel=L"\overline{T}_o")
plot!(plt_cop_AₒTₒ, xlabel=L"P(\overline{A}_o)", ylabel=L"P(\overline{T}_o)")

## C[F(tᶜ),F(Aₒ)]
JPD_tᶜAₒ, Cop_tᶜAₒ, plt_cop_tᶜAₒ, plt_cop_tᶜAₒ_OG = copula2d_fit_eval(GaussianCopula,alltᶜ,allAₒ,tᶜ_Fit,Aₒ_Fit,2^16,100)
### Add labels to the copula density plots
plot!(plt_cop_tᶜAₒ_OG, xlabel=L"\overline{t}^c", ylabel=L"\overline{A}_o")
plot!(plt_cop_tᶜAₒ, xlabel=L"P(\overline{t}^c)", ylabel=L"P(\overline{A}_o)")

### Sample from resulting joint probability distribution P(x₁,x₂)
# Binning of of the conditioning variable x₁
## Range for x₁ (lower and upper bound)   
x₁_hrange = mean(alltᶜ)+std(alltᶜ)   

Δtᶜ_μ = Δtᶜ_dist_params[1]
Δtᶜ_σ = Δtᶜ_dist_params[2]
Δtᶜ_var = (exp(Δtᶜ_σ^2)-1)*exp(2*Δtᶜ_μ+Δtᶜ_σ^2)
σᴮ = sqrt(Δtᶜ_var)    # Use of standard deviation of appropriate dataset to determine width of bin
wᴮ = 2*σᴮ   # Width of bin

x₁ˡᵇ = -(Int(floor(x₁_hrange/wᴮ))+0.5)*wᴮ   
x₁ᵘᵇ = (Int(floor(x₁_hrange/wᴮ))+0.5)*wᴮ 
Nᴮ = Int(floor((x₁ᵘᵇ-x₁ˡᵇ)/wᴮ)) + 1  # Number of bins
x₁_bins = range(x₁ˡᵇ,x₁ᵘᵇ,Nᴮ)   # Central values of bins for the x₁ range

plt_tcbins = plot(xlab=L"t~[s]", ylab=L"A_o~[m]")
for i ∈ 1:Nᴮ
    plot!(plt_tcbins, [x₁_bins[i].*Tₚⱼ; x₁_bins[i].*Tₚⱼ], [0; maximum(allAₒ)*Hₛⱼ], ls=:dash, leg=:false)
end

x₁ᶜ = x₁_bins[Int(floor(Nᴮ/2))+1] # Which bin?
hist_x2, hist_x1x2, x₁_samp, x₂_samp, x₁_range = conditional_sampling(JPD_tᶜAₒ,10000,x₁ᶜ,wᴮ)

dataset = x₂_samp;   dist_type = Weibull;  Nbins = 2^9
pdf_values, fitted_dist, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
plot!(hst_plt,  title="PDF - P(A̅ₒ|tᶜ bin = $(round(x₁ᶜ*1000)/1000)", xlab=L"\overline{A}_o", ylab=L"P(\overline{A}_o)")
FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])

Aₒ_WeiFits = zeros(Float64,Nᴮ,2)
Aₒ_WeiModes = zeros(Float64,Nᴮ)
plt_BinDist = plot(xlab=L"t^c", ylab=L"\overline{A}_o", zlab=L"P(\overline{A}_o)")
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
    plot3d!(plt_BinDist, x₁ᶜ*ones(Nbins), Aₒ_Wrang, pdf_values, lw=2, lab="t= $(round(x₁ᶜ*1000)/1000)")
    # plot3d!(hst_plt,  title="PDF - P(A̅ₒ|tᶜ bin = $(round(x₁ᶜ*1000)/1000)", xlab=L"\overline{A}_o", ylab=L"P(\overline{A}_o)")
    # display(hst_plt)

    Aₒ_WeiFits[j,:] = dist_params[:]
    kʷ = dist_params[1]
    λʷ = dist_params[2]
    mode_value = λʷ*((kʷ-1)/kʷ)^(1/kʷ)
    Aₒ_WeiModes[j] = mode_value
    # println("Weibull Mode Value: $mode_value")
end
display(plt_BinDist)

## 3-variate copula: C[F(tᶜ), F(Tₒ), F(Aₒ)]
JPD_tᶜTₒAₒ, Cop_tᶜTₒAₒ, plt_3vc_tᶜTₒ, plt_3vc_tᶜTₒ_OG, plt_3vc_tᶜAₒ, plt_3vc_tᶜAₒ_OG, plt_3vc_TₒAₒ, plt_3vc_TₒAₒ_OG = copula3d_fit_eval(GaussianCopula,alltᶜ,allTₒ, allAₒ,tᶜ_Fit,Tₒ_Fit, Aₒ_Fit,2^10,100)
### Add labels to the copula density plots
plot!(plt_3vc_tᶜTₒ_OG, xlabel=L"\overline{t}^c", ylabel=L"\overline{T}_o")
plot!(plt_3vc_tᶜAₒ_OG, xlabel=L"\overline{t}^c", ylabel=L"\overline{A}_o")
plot!(plt_3vc_TₒAₒ_OG, xlabel=L"\overline{T}_o", ylabel=L"\overline{A}_o")
plot!(plt_3vc_tᶜTₒ, xlabel=L"P(\overline{t}^c)", ylabel=L"P(\overline{T}_o)")
plot!(plt_3vc_tᶜAₒ, xlabel=L"P(\overline{t}^c)", ylabel=L"P(\overline{A}_o)")
plot!(plt_3vc_TₒAₒ, xlabel=L"P(\overline{T}_o)", ylabel=L"P(\overline{A}_o)")

plt_cop_tᶜTₒAₒ_OG = plot(plt_3vc_tᶜTₒ_OG, plt_3vc_tᶜAₒ_OG, plt_3vc_TₒAₒ_OG, layout= @layout[a b; c d], title = "3-variate copula f(tᶜ,T,A)")
plt_cop_tᶜTₒAₒ = plot(plt_3vc_tᶜTₒ, plt_3vc_tᶜAₒ, plt_3vc_TₒAₒ, layout= @layout[a b; c d], title = "3-variate copula f(tᶜ,T,A)")

#############################################################################################
# Heatmaps
### Tₒ-Aₒ
hplt_TA = histogram2d(allTₒ, allAₒ, normalize=:pdf, bins=:sqrt, show_empty_bins=true, color=:plasma, xlab = L"\overline{T}", ylab = L"\overline{A}")

### tᶜ-Aₒ
tcbins = range(-t̅ₑ, t̅ₑ, length=154)
Aobins = range(0, 1.2, length=25)
Tobins = range(0, 3.5, length=35)
hplt_tcA = histogram2d(alltᶜ, allAₒ, normalize=:pdf, bins=(tcbins,Aobins), show_empty_bins=true, color=:plasma, xlab = L"\overline{t}", ylab = L"\overline{A}")

hist = fit(Histogram, (alltᶜ, allAₒ), (tcbins,Aobins), closed=:right)
etc = hist.edges[1]
eA = hist.edges[2]
probabs = (hist.weights / sum(hist.weights))'
heat_tcA = heatmap(etc[1:end-1], eA[1:end-1], probabs, xlim=(-t̅ₑ,t̅ₑ), xlab = L"\overline{t}", ylab = L"\overline{A}")
# surface(etc[1:end-1], eA[1:end-1], probabs, xlim=(-t̅ₑ,t̅ₑ), xlab = L"\overline{t}", ylab = L"\overline{A}",color=:plasma)

### tᶜ-Tₒ
hplt_tcT = histogram2d(alltᶜ, allTₒ, normalize=:pdf, bins=(tcbins,Tobins), show_empty_bins=true, color=:plasma, xlab = L"\overline{t}", ylab = L"\overline{T}")

#############################################################################################
# Trying the tᶜ-A*tᶜ distribution
σᵗᶜ = std(alltᶜ)
alltᶜ_nz, iz = [], []
for i ∈ 1:length(alltᶜ)
    if alltᶜ[i] ≠ 0 && abs(alltᶜ[i]) < 3*σᵗᶜ
        push!(alltᶜ_nz, alltᶜ[i])
        push!(iz, i)
    end
end

alltᶜ_nz =  Float64.(alltᶜ_nz); iz = Int.(iz)
Aotc = allAₒ[iz]./alltᶜ_nz

Aotc_bins = range(minimum(Aotc), maximum(Aotc), length=Int(round(sqrt(length(Aotc)))))
hist_Aotc = fit(Histogram, Aotc, Aotc_bins, closed=:right)
p_Aotc = (hist_Aotc.weights / sum(hist_Aotc.weights))
plot(Aotc_bins[1:end-1], p_Aotc)

μᴬᵗ = mean(abs.(Aotc)); σᴬᵗ = std(abs.(Aotc))

Aotc_fin, ifin = [], []
for i ∈ 1:length(Aotc)
    if abs(Aotc[i]) < (Aotc_bins[2]-Aotc_bins[1])/2
        push!(Aotc_fin, Aotc[i])
        push!(ifin, i)
    end
end

Aotc_fin = Float64.(Aotc_fin)
alltᶜ_fin = alltᶜ_nz[ifin]

plot(alltᶜ_fin,Aotc_fin, seriestype=:scatter)#, ylim=(-1,1))
plot!(sort(alltᶜ_fin), π/2*exp.(-abs.(sort(alltᶜ_fin))).*coth.(sort(alltᶜ_fin)), lw=2)
# plot!(sort(alltᶜ_fin), exp(1) * sign.(sort(alltᶜ_fin)) .* exp.(-abs.(sort(alltᶜ_fin))), lw=2)
# plot!(sort(alltᶜ_fin), sign.(sort(alltᶜ_fin)) .* exp.(-sqrt.(abs.(sort(alltᶜ_fin)))), lw=2)

plot(alltᶜ,allAₒ, seriestype=:scatter)
plot!(alltᶜ_fin,Aotc_fin.*alltᶜ_fin, seriestype=:scatter)

histogram2d(alltᶜ_fin,Aotc_fin, normalize=:pdf, show_empty_bins=true, color=:plasma, bins=:sqrt)
histogram(Aotc_fin, normalize=:pdf, show_empty_bins=true, color=:plasma, bins=:sqrt)

dataset = Aotc_fin;   dist_type = Cauchy;  Nbins = 2^9
pdf_values, Aotc_Fit_fin, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
plot!(hst_plt,  title="PDF - P(A̅ₒ/t̅ᶜ)",  xlab=L"\frac{\overline{A}_o}{\overline{t}_c}", ylab=L"P(\frac{\overline{A}_o}{\overline{t}_c})")
FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])

dataset = alltᶜ_fin;   dist_type = Normal;  Nbins = 2^9
pdf_values, tᶜ_Fit_fin, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
plot!(hst_plt,  title="PDF - P(t̅ᶜ)", xlab=L"\overline{t}^c", ylab=L"P(\overline{t}^c)")
FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])

JPD_tᶜAtc, Cop_tᶜAtc, plt_cop_tᶜAtc, plt_cop_tᶜAtc_OG = copula2d_fit_eval(GaussianCopula,alltᶜ_fin,Aotc_fin,tᶜ_Fit_fin,Aotc_Fit_fin,2^10,100)
### Add labels to the copula density plots
plot!(plt_cop_tᶜAtc_OG, xlabel=L"\overline{t}^c", ylabel=L"\frac{\overline{A}_o}{\overline{t}_c}")

#############################################################################################
# Interpolation of most probable envelope shape from t-A heatmap
Aint = zeros(Float64, size(probabs)[2])

for i ∈ 1:size(probabs)[2]
    j = findmax(probabs[:,i])[2]
    Aint[i] = eA[j]
end

A_itp = interpolate(etc[1:end-1], Aint, BSplineOrder(4))
Aᵢ = A_itp.(t̅)

Aₒᵢ = zeros(Float64,Nₜ)
Aₒᵢ[Int(Nₜ/2-2^10):Int(Nₜ/2+2^10)] = Aᵢ[Int(Nₜ/2-2^10):Int(Nₜ/2+2^10)]

fitβ̃ = linear_fit(alltᶜ,allβ)
fitβ = linear_fit(alltᶜ.*T̃ₚᴿ,allβ)

# β̃reg = fitβ̃[1] .+ fitβ̃[2].*alltᶜ
β̃reg = fitβ̃[2].*alltᶜ
plot!(pltGEE5,alltᶜ,β̃reg)

#############################################################################################
# Gaussian Regression of mean envelope
Ø = Array{Float64}(undef,0)
Tₒ, Aₒ, tᶜₒ = [Ø[:] for _ = 1:3]
pltCONV = plot()

# Envelope Regression
if regress
    # Peak detection parameters
    MinPeakVal = 0 #1.25*std(η) 
    MinPeakDist = 2*Tcut
    # u = Hₛⱼ*Gmean
    u = Aₒᵢ*H̃ₛᴿ

    # Low-pass filtering of envelope
    # fcˡᵖ = 2*fcut;    fₛˡᵖ = round(1/dt)
    # Ntaps = nextpow(2,fₛˡᵖ/fcˡᵖ)   # No of taps
    # u = low_pass_filter(u,fcˡᵖ,fₛˡᵖ,Ntaps)

    # Runge-Kutta solution parameters (see in regress_gauss_fun.jl for details)
    ϵ, Nτ, dτ, Tlb, Tub = 1e-4, Int(1e4), 1e+0, 0*Tcut, Tₑ 

    # Call the Gaussian Regression function
    Tₒ, Aₒ, tᶜₒ, i⁺ₒ, G, pltCONV = regress_gauss_fun(t,u,MinPeakVal,MinPeakDist,ϵ,Nτ,dτ,Tlb,Tub)
else
    NoGEEs = 9
    u = Aₒᵢ*H̃ₛᴿ
    wdtc = rand(Δtᶜ_Fit,NoGEEs-1)
    tᶜₘ = 0.0
    # tᶜₘ = t̅[findmax(u)[2]]
    # # tᶜₘ = -rand()
    # tᶜₒ = [tᶜₘ-wdtc[1]; tᶜₘ; tᶜₘ+wdtc[2]].*T̃ₚᴿ
    tᶜₒ = [tᶜₘ-sum(wdtc[1:4]); tᶜₘ-sum(wdtc[1:3]); tᶜₘ-sum(wdtc[1:2]); tᶜₘ-wdtc[1]; tᶜₘ; tᶜₘ+wdtc[5]; tᶜₘ+sum(wdtc[5:6]); tᶜₘ+sum(wdtc[5:7]); tᶜₘ+sum(wdtc[5:8])].*T̃ₚᴿ
    
    # Conditional sampling for Aₒ
    Aₒ_range = range(0, maximum(allAₒ), length=2^9)
    plt_weibs = plot3d(xlab=L"t^c", ylab=L"A_o" ,zlab="Density")
    for i ∈ 1:NoGEEs
        # _, _, _, Aₒ_samp, _ = conditional_sampling(JPD_tᶜAₒ,10000,tᶜₒ[i]/T̃ₚᴿ,wᴮ)
        # pdf_Aₒ, fitted_dist, _ = dist_fit(Aₒ_samp, Weibull, 2^9)
        bin_tᶜ = Int(floor((tᶜₒ[i]/T̃ₚᴿ-x₁ˡᵇ)/wᴮ)+1)
        fitted_dist = Weibull(Aₒ_WeiFits[bin_tᶜ,1],Aₒ_WeiFits[bin_tᶜ,2])
        pdf_Aₒ = pdf.(fitted_dist, Aₒ_range)
        push!(Aₒ, rand(fitted_dist))
        plot3d!(plt_weibs, ones(2^9)*tᶜₒ[i], Aₒ_range,pdf_Aₒ, lw=2, lab=L"P(A_o)"*"for GEE $i")
    end
    display(plt_weibs)

    # Conditional sampling for Tₒ
    Tₒ_range = range(0, maximum(allTₒ), length=2^9)
    plt_weibs2 = plot3d(xlab=L"A_o", ylab=L"T_o" ,zlab="Density")
    for i ∈ 1:NoGEEs
        _, _, _, Tₒ_samp, _ = conditional_sampling(JPD_AₒTₒ,10000,Aₒ[i],std(allAₒ))
        pdf_Tₒ, fitted_dist, _ = dist_fit(Tₒ_samp, LogNormal, 2^9)
        Tₒⁱ = rand(fitted_dist)
        while Tₒⁱ < Tcut/T̃ₚᴿ
            Tₒⁱ = rand(fitted_dist)
        end
        push!(Tₒ, Tₒⁱ)
        plot3d!(plt_weibs2, ones(2^9)*Aₒ[i], Tₒ_range,pdf_Tₒ, lw=2, lab=L"P(T_o)"*"for GEE $i")
    end
    display(plt_weibs2)

    Aₒ = Aₒ.*H̃ₛᴿ
    Tₒ = Tₒ.*T̃ₚᴿ

    # i⁺ = Int.(round.(tᶜₒ./(t[2]-t[1]) .+ Int(Nₜ/2)))
    # Aₒ = u[i⁺]
    # Tₒ = rand(Tₒ_Fit,NoGEEs).*T̃ₚᴿ

    # Aₒ = sqrt(2)*exp.(-abs.(tᶜₒ./T̃ₚᴿ)).*coth.(tᶜₒ./T̃ₚᴿ) .* tᶜₒ./T̃ₚᴿ .*H̃ₛᴿ

    # P1 = rand(JPD_tᶜTₒAₒ,9)
    # tᶜₒ = P1[1,:].*T̃ₚᴿ
    # Tₒ = P1[2,:].*T̃ₚᴿ
    # Aₒ = P1[3,:].*H̃ₛᴿ

    # P2 = rand(JPD_TₒAₒ,9)
    # Tₒ = P2[1,:].*T̃ₚᴿ
    # Aₒ = P2[2,:].*H̃ₛᴿ

    G = gauss_fun(t, Aₒ, tᶜₒ, Tₒ)
end

plot!(plt_cop_TₒAₒ_OG, [Tₒ./T̃ₚᴿ], [Aₒ./H̃ₛᴿ], seriestype=:scatter, lab="Sampled")
plot!(plt_cop_tᶜAₒ_OG, [tᶜₒ./T̃ₚᴿ], [Aₒ./H̃ₛᴿ], seriestype=:scatter, lab="Sampled")
plot!(plt_tcbins, tᶜₒ, Aₒ, xlim=(-64,64), seriestype=:scatter, lab="Sampled")
display(plt_tcbins)

# Reconstruction
lenT = length(Tₒ)
Nₛ₂ = Int(nextpow(2,Nₜ)/2)
TOT = zeros(ComplexF64, Nₜ)  # Total surface elevation
SP_TOT = zeros(Float64, Nₛ₂) # Total spectrum of propagated Gaussian EWGs

plt_gn = plot(palette=:darkrainbow)
pltEWG = plot(palette=:darkrainbow)
pltSPEC = plot(palette=:darkrainbow)
plt_demo = plot(palette=[cb[8];cb[11];cb[4]])

prop = 0
β̃ = fitβ[2]*tᶜₒ
Ωₙ = rand(ΩFit,lenT)
ωₚ = 2π/Tₚⱼ
ωᵧ = 4*acos(exp(-1/4))./Tₒ
ω₁ = mean(allΩ)*ωₚ
ω₂ = 0.0
for n ∈ 1:lenT
    # Propagated EWG
    global ω₂ = ω₁ - ωᵧ[n]
    # gₙ, ηᶜ, ηₙ, FR, MAG = recon(prop, Aₒ[n], tᶜₒ[n], Tₒ[n], β̃[n], 1, 1, mean(allΩ)/Tₚⱼ, t)
    # gₙ, ηᶜ, ηₙ, FR, MAG = recon(prop, Aₒ[n]/Hₛⱼ, tᶜₒ[n]/Tₚⱼ, Tₒ[n]/Tₚⱼ, β̃[n], Hₛⱼ, Tₚⱼ, mean(allΩ), t)
    gₙ, ηᶜ, ηₙ, FR, MAG = recon(prop, Aₒ[n]/H̃ₛᴿ, tᶜₒ[n]/T̃ₚᴿ, Tₒ[n]/T̃ₚᴿ, β̃[n], H̃ₛᴿ, T̃ₚᴿ, Ωₙ[n], t)

    TOT[:] = TOT[:] .+ real.(ηₙ)
    SP_TOT[:] = SP_TOT[:] .+ MAG

    plot!(plt_gn, t, gₙ, xlab = "t [s]", ylab = "A [m]", lab = "g$(n)", lw=2, legend=:outerbottom, legendcolumns=min(lenT,8))
    plot!(pltEWG, t, real(ηₙ), xlab = "t [s]", ylab = "η [m]", lab = "η$(n)", lw=2, legend=:outerbottom, legendcolumns=min(lenT,8))
    plot!(pltSPEC, FR, MAG, xlab = L"f [Hz]", ylab = L"A [m]", lab = "g$(n)", lw=2, legend=:outerbottom, legendcolumns=min(lenT,8))
    plot!(pltSPEC, xlim=(0,fcut))

    if n == 1
        plot!(plt_demo, t, gₙ, xlab = "t [s]", ylab = "[m]", lab = L"g_1", lw=2)
        plot!(plt_demo, t, real(ηᶜ), lab="Temporal term", lw=2)
        plot!(plt_demo, t, real(ηₙ), lab="Propagated elementary envelope", lw=2)
        # plot!(plt_demo, xlim=(-4*T̅ₒ[n]*Tₚ,4*T̅ₒ[n]*Tₚ))
    end                         
end
# _, _, η₂, _ = recon(prop, 0.1, 0, 3, -π, Hₛⱼ, Tₚⱼ, 2π/Tᵖ₂₋, t)
ηᴳ = real(TOT) #.+ real(η₂)
# ηᴳ = real(TOT) .+ η₂₋
plt_etaG = plot(t,ηᴳ,xlab="t [s]", ylab="η [m]", lw=2, xlim=(t[Int(Nₜ/2-2^10)], t[Int(Nₜ/2+2^10)]))

# Gaussian approximation against envelope and scatter plot of peaks
plt_envs = plot(xlab = L"t~[s]", ylab = L"[m]", palette=[cb[4];cb[8];cb[11]], legend=:outertop, legendcolumns=4)
if regress
    plot!(t, [Hₛⱼ*Gmean u G], lab = [L"u(t)" L"u_{flt}(t)" L"G(t)"], lw=[2 2 2])
else
    plot!(t, [Aᵢ*H̃ₛᴿ u], lab = [L"u(t)" L"u_{flt}(t)" L"G(t)"], lw=[2 2 2])
end
plot!(tᶜₒ, Aₒ, seriestype=:scatter, ms=2, mc=:red, lab = "Peaks")

if wout
    fid = joinpath(homedir(),"OpenFAST_Simulations","ExtElev",case_id,"DWE.Elev")
    open(fid, "w")
    writedlm(fid, [t.+t[end] ηᴳ], '\t')
end

stats_path = joinpath(libpath,"SE",case_id,"0")

#############################################################################################
# PLOTS
for i ∈ 1:length(AllFitPlots)
    display(AllFitPlots[i])
    savefig(AllFitPlots[i],joinpath(stats_path,"$(plt_img_names[i]).png"))
    savefig(AllFitPlots[i],joinpath(stats_path,"$(plt_img_names[i]).svg"))
end

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

# display(pltSP2)
# savefig(pltSP2, joinpath(stats_path,"sp_2nd-.png"))
# savefig(pltSP2, joinpath(stats_path,"sp_2nd-.svg"))
# display(pltPHI2)
# savefig(pltPHI2, joinpath(stats_path,"phi_2nd-.png"))
# savefig(pltPHI2, joinpath(stats_path,"phi_2nd-.svg"))
# display(plt_eta2)
# savefig(plt_eta2, joinpath(stats_path,"eta_2nd-.png"))
# savefig(plt_eta2, joinpath(stats_path,"eta_2nd-.svg"))
display(hplt_TA)
savefig(hplt_TA, joinpath(stats_path,"hist_To-Ao.png"))
savefig(hplt_TA, joinpath(stats_path,"hist_To-Ao.svg"))
display(hplt_tcA)
savefig(hplt_tcA, joinpath(stats_path,"hist_tc-Ao.png"))
savefig(hplt_tcA, joinpath(stats_path,"hist_tc-Ao.svg"))
display(hplt_tcT)
savefig(hplt_tcT, joinpath(stats_path,"hist_tc-To.png"))
savefig(hplt_tcT, joinpath(stats_path,"hist_tc-To.svg"))
display(plt_cop_TₒAₒ_OG)
savefig(plt_cop_TₒAₒ_OG, joinpath(stats_path,"Ao-To_copula.png"))
savefig(plt_cop_TₒAₒ_OG, joinpath(stats_path,"Ao-To_copula.svg"))
display(plt_cop_tᶜAₒ_OG)
savefig(plt_cop_tᶜAₒ_OG, joinpath(stats_path,"tc-Ao_copula.png"))
savefig(plt_cop_tᶜAₒ_OG, joinpath(stats_path,"tc-Ao_copula.svg"))
display(plt_cop_tᶜTₒAₒ_OG)
savefig(plt_cop_tᶜTₒAₒ_OG, joinpath(stats_path,"var3_copula.png"))
savefig(plt_cop_tᶜTₒAₒ_OG, joinpath(stats_path,"var3_copula.svg"))
# display(pltG)
# savefig(pltG, joinpath(stats_path,"envelopes.png"))
# savefig(pltG, joinpath(stats_path,"envelopes.svg"))
display(plt_envs)
savefig(plt_envs, joinpath(stats_path,"env_reg.png"))
savefig(plt_envs, joinpath(stats_path,"env_reg.svg"))
display(plt_gn)
savefig(plt_gn, joinpath(stats_path,"Eenvs.png"))
savefig(plt_gn, joinpath(stats_path,"Eenvs.svg"))
display(pltEWG)
savefig(pltEWG, joinpath(stats_path,"EWGs.png"))
savefig(pltEWG, joinpath(stats_path,"EWGs.svg"))
display(pltSPEC)
savefig(pltSPEC, joinpath(stats_path,"EWGs_spec.png"))
savefig(pltSPEC, joinpath(stats_path,"EWGs_spec.svg"))
# display(plt_demo)
# savefig(plt_demo, joinpath(stats_path,"carrier.png"))
# savefig(plt_demo, joinpath(stats_path,"carrier.svg"))
display(plt_etaG)
savefig(plt_etaG, joinpath(stats_path,"DWE.png"))
savefig(plt_etaG, joinpath(stats_path,"DWE.svg"))

if regress
    display(pltCONV)
    savefig(pltCONV, joinpath(stats_path,"convergence.png"))
    savefig(pltCONV, joinpath(stats_path,"convergence.svg"))
end


# Summarise stats from all simulations
TotΩ = allΩ;  Totβ = allβ;  Totβ̇ = allβ̇;  TotN = allN
TotAₒ = allAₒ;    TotTₒ = allTₒ;    Tottᶜ = alltᶜ
TotΔtᶜ = allΔtᶜ;  TotΔTₑᵥ = allΔTₑᵥ

# TotΩ = [TotΩ; allΩ];  Totβ = [Totβ; allβ];  Totβ̇ = [Totβ̇; allβ̇];  TotN = [TotN;allN]
# TotAₒ = [TotAₒ; allAₒ];    TotTₒ = [TotTₒ; allTₒ];    Tottᶜ = [Tottᶜ; alltᶜ]
# TotΔtᶜ = [TotΔtᶜ; allΔtᶜ];  TotΔTₑᵥ = [TotΔTₑᵥ; allΔTₑᵥ] 

# # Heatmaps
# ### Tₒ-Aₒ
# hplt_TA = histogram2d(TotTₒ, TotAₒ, normalize=:pdf, bins=:sqrt, show_empty_bins=true, color=:plasma, xlab = L"\overline{T}", ylab = L"\overline{A}")

### tᶜ-Aₒ
# tcbins = range(-t̅ₑ, t̅ₑ, length=74)
# Aobins = range(0, 1.2, length=35)
# Tobins = range(0, 3.5, length=35)
# hplt_tcA = histogram2d(Tottᶜ, TotAₒ, normalize=:pdf, bins=(tcbins,Aobins), show_empty_bins=true, color=:plasma, xlab = L"\overline{t}", ylab = L"\overline{A}")

# hist = fit(Histogram, (Tottᶜ, TotAₒ), (tcbins,Aobins), closed=:right)
# etc = hist.edges[1]
# eA = hist.edges[2]
# probabs = (hist.weights / sum(hist.weights))'
# heat_tcA = heatmap(etc[1:end-1], eA[1:end-1], probabs, xlim=(-t̅ₑ,t̅ₑ), xlab = L"\overline{t}", ylab = L"\overline{A}")

# # Statistical Fits
# dataset = TotΩ;   dist_type = LogNormal;  Nbins = 2^9
# pdf_values, ΩFit, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
# plot!(hst_plt, title="PDF - P(Ω)", xlab=L"Ω", ylab=L"P(Ω)")
# FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])

# ### β̇
# dataset = Totβ̇ ;   dist_type = Normal;  Nbins = 2^9
# pdf_values, fitted_dist, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
# plot!(hst_plt, title="PDF - P(β̇)", xlab=L"\dot{\beta}~[rad]", ylab=L"P(\dot{\beta})")
# FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])

# ### N
# dataset = TotN;   dist_type = Normal;  Nbins = 2^9
# pdf_values, fitted_dist, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
# plot!(hst_plt, title="PDF: P(N̅)", xlab=L"\overline{N}", ylab=L"P(\overline{N})")
# FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])

# ### Δtᶜ
# dataset = TotΔtᶜ;   dist_type = LogNormal;  Nbins = 2^9
# pdf_values, Δtᶜ_Fit, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
# plot!(hst_plt, title="PDF - P(Δtᶜ)", xlab=L"\overline{\Delta t^c}", ylab=L"P(\overline{\Delta t^c})")
# FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])

# ### Tₒ
# dataset = TotTₒ;   dist_type = LogNormal;  Nbins = 2^9
# pdf_values, Tₒ_Fit, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
# plot!(hst_plt, title="PDF - P(T̅ₒ)", xlab=L"\overline{T}_o", ylab=L"P(\overline{T}_o)")
# FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])

# ### tᶜ
# dataset = Tottᶜ;   dist_type = Cauchy;  Nbins = 2^9
# pdf_values, tᶜ_Fit, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
# plot!(hst_plt,  title="PDF - P(t̅ᶜ)", xlab=L"\overline{t}^c", ylab=L"P(\overline{t}^c)")
# FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])

# ### Aₒ
# dataset = TotAₒ;   dist_type = Weibull;  Nbins = 2^9
# pdf_values, Aₒ_Fit, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
# plot!(hst_plt,  title="PDF - P(A̅ₒ)", xlab=L"\overline{A}_o", ylab=L"P(\overline{A}_o)")
# FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])

# # Copulas
# ## C[F(Tₒ),F(Aₒ)]
# JPD_TₒAₒ, Cop_TₒAₒ, plt_cop_TₒAₒ, plt_cop_TₒAₒ_OG = copula2d_fit_eval(GaussianCopula,TotTₒ,TotAₒ,Tₒ_Fit,Aₒ_Fit,2^16,100)
# ### Add labels to the copula density plots
# plot!(plt_cop_TₒAₒ_OG, xlabel=L"\overline{T}_o", ylabel=L"\overline{A}_o")
# plot!(plt_cop_TₒAₒ, xlabel=L"P(\overline{T}_o)", ylabel=L"P(\overline{A}_o)")

# ## C[F(tᶜ),F(Aₒ)]
# JPD_tᶜAₒ, Cop_tᶜAₒ, plt_cop_tᶜAₒ, plt_cop_tᶜAₒ_OG = copula2d_fit_eval(GaussianCopula,Tottᶜ,TotAₒ,tᶜ_Fit,Aₒ_Fit,2^16,100)
# ### Add labels to the copula density plots
# plot!(plt_cop_tᶜAₒ_OG, xlabel=L"\overline{t}^c", ylabel=L"\overline{A}_o")
# plot!(plt_cop_tᶜAₒ, xlabel=L"P(\overline{t}^c)", ylabel=L"P(\overline{A}_o)")
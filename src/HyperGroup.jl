module HyperGroup
# Load necessary modules
using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles
using FFTW, DSP, Statistics, BSplineKit, StatsBase, Distributions, Copulas, Optim
using CurveFit: linear_fit
using StatsPlots
import ColorSchemes.darkrainbow
import Polynomials: fit as pfit

gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))
cb = darkrainbow

include("WAVEGEM.jl")
import .WAVEGEM

# Include necessary scripts for functions
include("func/text_process.jl"); include("func/signal_processing.jl")
include("func/peak_detect.jl");  include("func/runge_kutta.jl")
include("func/fugees.jl");       include("func/directories.jl")
include("func/wave_theory.jl");  include("func/jonswap.jl")
include("func/stat_process.jl")

#############################################################################################
# Global Variables
## Paths
case_id, _, rundir = WAVEGEM.GlobPaths[2:4]
OFASTpath = WAVEGEM.GlobPaths[end]
## Other
stats_path, gen_method, wout = WAVEGEM.HyperGroupInp[1:3]
ρ, g, _, _, Hₛⱼ, Tₚⱼ, γ, fᶜᵘᵗ = WAVEGEM.GlobInp0
Tₑ = WAVEGEM.Tₑ
T₀ₛ, T₀ₕ, T₀ₚ = WAVEGEM.FOWT[:]

Tᶜᵘᵗ = 1/fᶜᵘᵗ   # Cut-off periodv
ϵₛₛ = Hₛⱼ/Tₚⱼ^2 * 2π/g  # Effective steepness

# Some initialisations
Ø = Array{Float64}(undef,0)
Øᵃⁿʸ = Array{Any}(undef,0)
AllPlots = Øᵃⁿʸ[:]
plt_img_names = Array{String}(undef,0)

#############################################################################################
# Load GEE datasets
fid = joinpath(stats_path,"summary_GEE_datasets")
GEE_dsets = parse_fxw(fid, 1) 
## Assign
allA̅, allT̅, allT̅ᵣ, allT̅ₑ, allt̅ᶜ, allt̅ᶜᵣ,  allt̅ᶜₑ, allβ = [GEE_dsets[:,i] for i ∈ 1:8]
maxA̅ = maximum(allA̅);   maxT̅ᵣ = maximum(allT̅ᵣ)

# Load Δt datasets
fid = joinpath(stats_path,"summary_Δt_datasets")
Δt_dsets = parse_fxw(fid, 1) 
## Assign
allΔt̅ᶜ, allΔt̅ᶜᵣ, allΔt̅ᶜₑ = [Δt_dsets[:,i] for i ∈ 1:3]

# Read dataset statistics (x̅ s Med	IQR	Skew Kurt)
fid = joinpath(stats_path,"sample_stats")
dsets_stats = parse_txt_flex(fid, 1, 1) 
## Assign 
A̅_dstat,T̅_dstat,T̅ᵣ_dstat, T̅ₑ_dstat, t̅ᶜ_dstat, t̅ᶜᵣ_dstat, t̅ᶜₑ_dstat = [dsets_stats[i,:] for i ∈ 1:7]  # GEE dataset statistics

# Read Distribution parameters
fid = joinpath(stats_path,"Dist_stats")
dist_stats = parse_txt_flex(fid, 1, 2) 
## Assign
A̅_dist,T̅_dist,T̅ᵣ_dist, T̅ₑ_dist, t̅ᶜ_dist, t̅ᶜᵣ_dist, t̅ᶜₑ_dist = [dist_stats[i,:] for i ∈ 1:7]  # GEE distributions
Ω_dist, β̇_dist, ΔTₑᵥ_dist, N̅_dist = [dist_stats[i,:] for i ∈ 9:12]  # CWE distributions
Δtᶜᵣ_dist = dist_stats[end-1,:]
Δtᶜₑ_dist = dist_stats[end,:]

## Define fiitted distributions
A̅_Fit = Weibull(A̅_dist[1],A̅_dist[2])
T̅ᵣ_Fit = LogNormal(T̅ᵣ_dist[1],T̅ᵣ_dist[2])
T̅ₑ_Fit = LogNormal(T̅ₑ_dist[1],T̅ₑ_dist[2])
t̅ᶜ_Fit = Cauchy(t̅ᶜ_dist[1],t̅ᶜ_dist[2])
t̅ᶜᵣ_Fit = Cauchy(t̅ᶜᵣ_dist[1],t̅ᶜᵣ_dist[2])
t̅ᶜₑ_Fit = Normal(t̅ᶜₑ_dist[1],t̅ᶜₑ_dist[2])
Ω_Fit = LogNormal(Ω_dist[1],Ω_dist[2])
β̇_Fit = Normal(β̇_dist[1],β̇_dist[2])
ΔTₑᵥ_Fit = LogNormal(ΔTₑᵥ_dist[1],ΔTₑᵥ_dist[2])
N̅_Fit = LogNormal(N̅_dist[1],N̅_dist[2])
Δtᶜᵣ_Fit = LogNormal(Δtᶜᵣ_dist[1],Δtᶜᵣ_dist[2])
Δtᶜₑ_Fit = LogNormal(Δtᶜₑ_dist[1],Δtᶜₑ_dist[2])

fitβ = linear_fit(allt̅ᶜᵣ.*T₀ₚ,allβ)
# plot(allt̅ᶜᵣ.*T₀ₚ, allβ, seriestype=:scatter)
# plot!(allt̅ᶜᵣ.*T₀ₚ, fitβ[2].*allt̅ᶜᵣ.*T₀ₚ .+ fitβ[1])

# beta = mod2pi.(fitβ[2].*t .+ fitβ[1])
# plot(t,beta, xlim=(-12,12))
#############################################################################################
# Determine number of Gaussian Elementary Envelopes
ΔTₑᵥᴴ = nextpow(2,ΔTₑᵥ_dist[5]*(ϵₛₛ*T₀ₛ)) / (ϵₛₛ*T₀ₛ)
NoGEEs = Int(floor(N̅_dist[5] * ΔTₑᵥᴴ * ϵₛₛ) + 1)
if mod(NoGEEs,2) == 0
    NoGEEs -= 1
end

# Time vector (non-dimensional)
pow2 = log(nextpow(2,ΔTₑᵥ_dist[5]*(ϵₛₛ*T₀ₛ)))/log(2) + 1
t̅ₑ = 2^pow2/T₀ₚ        # Duration
Nₜ::Int64 = 2^12    # No of time steps

dt = t̅ₑ/(Nₜ-1)  # Time step
t̅ = zeros(Float64,Nₜ)   # Initialise non-dimensional time vector
[t̅[i] = (i-1)*dt-t̅ₑ/2 for i ∈ 1:Nₜ] # Non-dimensional time vector
t = t̅*T₀ₚ   # Dimensional time vector

#############################################################################################
# Some heatmaps
hplt_tᶜᵣA, P₁₂, e₁, e₂, W¹² = joint_pdf_hmap(allt̅ᶜᵣ, allA̅, 1)
plot!(hplt_tᶜᵣA, xlab = L"\overline{t}_d", ylab = L"\overline{A}", title=L"P(\overline{t}_d,\overline{A})", xlim=(-t̅ₑ,t̅ₑ), ylim=(minimum(allA̅),1))
# push!(AllPlots,hplt_tᶜᵣA);    push!(plt_img_names,"hmap_trA")

hplt_ATᵣ, P₂₃, e₂, e₃,_ = joint_pdf_hmap(allA̅, allT̅ᵣ, 1)
hmap_PTAt = contourf(e₁[1:end-1],e₃[1:end-1],(P₂₃*100)*(P₁₂*100), show_empty_bins=true, color=:plasma, colorbar_title="P[%]",lw=0)
plot!(hmap_PTAt, title=L"P(T|t) = E_A[P(T|A)|t]", xlab=L"\overline{t}_d", ylab=L"\overline{T}_d",xlim=(-t̅ₑ,t̅ₑ), ylim=(minimum(allT̅ᵣ),1))
push!(AllPlots,hmap_PTAt);    push!(plt_img_names,"TAt_conditional")

hplt_TᵣA,_,_,_,W³² = joint_pdf_hmap(allT̅ᵣ, allA̅, 1)
plot!(xlab = L"\overline{T}_d", ylab = L"\overline{A}", title=L"P(\overline{T}_d,\overline{A})", xlim=(minimum(allT̅ᵣ),1), ylim=(minimum(allA̅),1))

hplt_tᶜₑA, = joint_pdf_hmap(allt̅ᶜₑ, allA̅, 1) 
plot!(hplt_tᶜₑA,xlab = L"\overline{t}_{ev}", ylab = L"\overline{A}", title=L"P(\overline{t}_{ev},\overline{A})", xlim=(-1,1), ylim=(minimum(allA̅),1))

hplt_tᶜᵣTᵣ,P₁₃,_,_,W¹³ = joint_pdf_hmap(allt̅ᶜᵣ, allT̅ᵣ, 1)
plot!(hplt_tᶜᵣTᵣ, xlab = L"\overline{t}_d", ylab = L"\overline{T}_d", title=L"P(\overline{t}_d,\overline{T}_d)", xlim=(-t̅ₑ,t̅ₑ), ylim=(minimum(allT̅ᵣ),1))
# push!(AllPlots,hplt_tᶜᵣTᵣ);    push!(plt_img_names,"hmap_P(td,Td)")

hplt_tᶜₑTᵣ,_ = joint_pdf_hmap(allt̅ᶜₑ, allT̅ᵣ, 1)
plot!(hplt_tᶜₑTᵣ, xlab = L"\overline{t}_{ev}", ylab = L"\overline{T}_d", title=L"P(\overline{t}_{ev},\overline{T}_d)", xlim=(-1,1), ylim=(minimum(allT̅ᵣ),1))

hmap_TᵣA, P₃₂, e₂, e₃,_ = joint_pdf_hmap(allT̅ᵣ, allA̅, 1)
hmap_PATt = contourf(e₁[1:end-1],e₃[1:end-1],(P₃₂*100)*(P₁₃*100), show_empty_bins=true, color=:plasma, colorbar_title="P[%]",lw=0)
plot!(hmap_PATt, title=L"P(A|t) =  E_T[P(A|T)|t]", xlab=L"\overline{t}_d", ylab=L"\overline{A}", xlim=(-t̅ₑ,t̅ₑ), ylim=(minimum(allA̅),1))
push!(AllPlots,hmap_PATt);    push!(plt_img_names,"ATt_coinditional")

#############################################################################################
## Copula: C[F(Tᵣ),F(A)]
JPD_TᵣA, cop_TᵣA, plt_cop_TᵣA, plt_cop_TᵣA_OG = copula2d_fit_eval(GaussianCopula,allT̅ᵣ,allA̅,T̅ᵣ_Fit,A̅_Fit,2^16,100)
### Add labels to the copula density plots
plot!(plt_cop_TᵣA_OG, xlabel=L"\overline{T}", ylabel=L"\overline{A}")
plot!(plt_cop_TᵣA, xlabel=L"P(\overline{T})", ylabel=L"P(\overline{A})")
cop_h2D_TᵣA,_ = copula_2D_hmap(JPD_TᵣA,2^16)
plot!(cop_h2D_TᵣA, xlab = L"\overline{T}_d", ylab = L"\overline{A}", title="Copula Density", xlim=(minimum(allT̅ᵣ),1), ylim=(minimum(allA̅),1))

## Copula: C[F(A),F(Tᵣ)]
JPD_ATᵣ, Cop_ATᵣ, plt_cop_ATᵣ, plt_cop_ATᵣ_OG = copula2d_fit_eval(GaussianCopula,allA̅,allT̅ᵣ,A̅_Fit,T̅ᵣ_Fit,2^16,100)
### Add labels to the copula density plots
plot!(plt_cop_ATᵣ_OG, xlabel=L"\overline{A}", ylabel=L"\overline{T}_d")
plot!(plt_cop_ATᵣ, xlabel=L"P(\overline{A})", ylabel=L"P(\overline{T}_d)")

## Copula: C[F(tᶜᵣ),F(Aₒ)]
JPD_tᶜᵣA, Cop_tᶜᵣA, plt_cop_tᶜᵣA, plt_cop_tᶜᵣA_OG = copula2d_fit_eval(GaussianCopula,allt̅ᶜᵣ,allA̅,t̅ᶜᵣ_Fit,A̅_Fit,2^16,100)
### Add labels to the copula density plots
plot!(plt_cop_tᶜᵣA_OG, xlabel=L"\overline{t}_d", ylabel=L"\overline{A}")
plot!(plt_cop_tᶜᵣA, xlabel=L"P(\overline{t}_d)", ylabel=L"P(\overline{A})")
_,x₁_samp, x₂_samp = copula_2D_histogram(JPD_tᶜᵣA,2^16)
i_out = findall(abs.(x₁_samp) .< 5*std(allt̅ᶜᵣ));    tᶜᵣ_samp = x₁_samp[i_out];  A_samp = x₂_samp[i_out]
cop_h2D_tᶜᵣA,_ = joint_pdf_hmap(tᶜᵣ_samp, A_samp, 1)
plot!(cop_h2D_tᶜᵣA, xlab = L"\overline{t}_d", ylab = L"\overline{A}", title="Copula Density", xlim=(-t̅ₑ,t̅ₑ), ylim=(0,1))

## Copula: C[F(tᶜₑ),F(A)]
JPD_tᶜₑA, Cop_tᶜₑA, plt_cop_tᶜₑA, plt_cop_tᶜₑA_OG = copula2d_fit_eval(GaussianCopula,allt̅ᶜₑ,allA̅,t̅ᶜₑ_Fit,A̅_Fit,2^16,100)
### Add labels to the copula density plots
plot!(plt_cop_tᶜₑA_OG, xlabel=L"\overline{t}_{ev}", ylabel=L"\overline{A}")
plot!(plt_cop_tᶜₑA, xlabel=L"P(\overline{t}_{ev})", ylabel=L"P(\overline{A})")
cop_h2D_tᶜₑA,_ = copula_2D_hmap(JPD_tᶜₑA,2^16)
plot!(cop_h2D_tᶜₑA, xlab = L"\overline{t}_{ev}", ylab = L"\overline{A}", title="Copula Density", xlim=(-1,1), ylim=(0,1))

## Copula: C[F(tᶜᵣ),F(Tᵣ)]
JPD_tᶜᵣTᵣ, Cop_tᶜᵣTᵣ, plt_cop_tᶜᵣTᵣ, plt_cop_tᶜᵣTᵣ_OG = copula2d_fit_eval(GaussianCopula,allt̅ᶜᵣ,allT̅ᵣ,t̅ᶜᵣ_Fit,T̅ᵣ_Fit,2^16,100)
### Add labels to the copula density plots
plot!(plt_cop_tᶜᵣTᵣ_OG, xlabel=L"\overline{t}_d", ylabel=L"\overline{T}_d")
plot!(plt_cop_tᶜᵣTᵣ, xlabel=L"P(\overline{t}_d)", ylabel=L"P(\overline{T}_d)")
_,x₁_samp, x₂_samp = copula_2D_histogram(JPD_tᶜᵣTᵣ,2^16)
i_out = findall(abs.(x₁_samp) .< 5*std(allt̅ᶜᵣ));    tᶜᵣ_samp = x₁_samp[i_out];  Tᵣ_samp = x₂_samp[i_out]
cop_h2D_tᶜᵣTᵣ,_ = joint_pdf_hmap(tᶜᵣ_samp, Tᵣ_samp, 1)
plot!(cop_h2D_tᶜᵣTᵣ, xlab = L"\overline{t}_d", ylab = L"\overline{T}_d", title="Copula Density", xlim=(-t̅ₑ,t̅ₑ), ylim=(minimum(allT̅ᵣ),1))

## Copula: C[F(tᶜₑ),F(Tᵣ)]
JPD_tᶜₑTᵣ, Cop_tᶜₑTᵣ, plt_cop_tᶜₑTᵣ, plt_cop_tᶜₑTᵣ_OG = copula2d_fit_eval(GaussianCopula,allt̅ᶜₑ,allT̅ᵣ,t̅ᶜₑ_Fit,T̅ᵣ_Fit,2^16,100)
### Add labels to the copula density plots
plot!(plt_cop_tᶜₑTᵣ_OG, xlabel=L"\overline{t}_{ev}", ylabel=L"\overline{T}_d")
plot!(plt_cop_tᶜₑTᵣ, xlabel=L"P(\overline{t}_{ev})", ylabel=L"P(\overline{T}_d)")
cop_h2D_tᶜₑTᵣ,_ = copula_2D_hmap(JPD_tᶜₑTᵣ,2^16)
plot!(cop_h2D_tᶜₑTᵣ, xlab = L"\overline{t}_{ev}", ylab = L"\overline{T}_d", title="Copula Density", xlim=(-1,1), ylim=(minimum(allT̅ᵣ),1))

mplt_tᶜᵣTᵣA = plot(hplt_TᵣA, cop_h2D_TᵣA, hplt_tᶜᵣA, cop_h2D_tᶜᵣA, layout=@layout[a b; c d])
push!(AllPlots,mplt_tᶜᵣTᵣA);    push!(plt_img_names,"multi_hmap1")

mplt_tT = plot(hplt_tᶜᵣTᵣ, cop_h2D_tᶜᵣTᵣ, hplt_tᶜₑTᵣ, cop_h2D_tᶜₑTᵣ, layout=@layout[a b; c d])
push!(AllPlots,mplt_tT);    push!(plt_img_names,"multi_hmap2")

mplt_tA = plot(hplt_tᶜᵣA, cop_h2D_tᶜᵣA, hplt_tᶜₑA, cop_h2D_tᶜₑA, layout=@layout[a b; c d])
push!(AllPlots,mplt_tA);    push!(plt_img_names,"multi_hmap3")

#############################################################################################
# Distributions of A̅ by binning joint probability distribution P(T̅ᵣ,A̅)
T̅ᵣ_bins, A̅_WeiFits, A̅_WeiModes, plt_Tr_bins, plt_BinDist0 = distros_from_jpd(allT̅ᵣ, allA̅, Weibull,L"\overline{T}_d")
plot!(plt_BinDist0, xlab=L"\overline{A}", ylab=L"P(\overline{A})")
push!(AllPlots,plt_BinDist0);    push!(plt_img_names,"ATr_BinDists")

# Interpolation of shape and scale parameters of resulting Weibull distributions per bin
## B-spline interpolation for κ
κ_itp = interpolate(T̅ᵣ_bins, A̅_WeiFits[:,1], BSplineOrder(4))    
T̅ᵣ_range = range(T̅ᵣ_bins[1],T̅ᵣ_bins[end],1000)
κⁱ = κ_itp.(T̅ᵣ_range)  

plt_kLN = plot(title="A̅ Distributions: Weibull shape parameter", xlab=L"\overline{T}_d", ylab=L"\kappa", palette=[cb[11];cb[8];cb[4]])
plot!(plt_kLN, T̅ᵣ_bins, A̅_WeiFits[:,1], seriestype=:scatter, lab=L"\kappa"*"(bins)")
plot!(plt_kLN, T̅ᵣ_range, κⁱ, lab=L"\kappa (\overline{T}_d)"*" B-Spline interpolation", lw=2)

## Polynomial interpolation
k_poly = pfit(T̅ᵣ_bins, A̅_WeiFits[:,1],4); k_coefs = round.(k_poly[:]*1000)./1000
k_eval = k_poly.(T̅ᵣ_range)

plot!(plt_kLN, T̅ᵣ_range, k_eval, lw=2)#, lab=L"\kappa (\overline{T}_d)"*"=$(k_coefs[5])"*L"x^4"*"$(k_coefs[4])"*L"x^3"*"+$(k_coefs[3])"*L"x^2"*"$(k_coefs[2])"*L"x"*"$(k_coefs[1])")
push!(AllPlots,plt_kLN);    push!(plt_img_names,"k_ATr_poly")

## B-spline interpolation for λ
λ_itp = interpolate(T̅ᵣ_bins, A̅_WeiFits[:,2], BSplineOrder(4))    
λⁱ = λ_itp.(T̅ᵣ_range)  

plt_λLN = plot(title="A̅ Distributions: Weibull scale parameter", xlab=L"\overline{T}_d", ylab=L"\lambda", palette=[cb[11];cb[8];cb[4]])
plot!(plt_λLN, T̅ᵣ_bins, A̅_WeiFits[:,2], seriestype=:scatter, lab=L"\lambda"*"(bins)")
plot!(plt_λLN, T̅ᵣ_range, λⁱ, lab=L"\lambda (\overline{T}_d)"*" B-Spline interpolation", lw=2)

## Polynomial interpolation
λ_poly = pfit(T̅ᵣ_bins, A̅_WeiFits[:,2],3); λ_coefs = round.(λ_poly[:]*1000)./1000
λ_eval = λ_poly.(T̅ᵣ_range)

plot!(plt_λLN, T̅ᵣ_range, λ_eval, lw=2)#, lab=L"\lambda (\overline{T})"*"=$(λ_coefs[5])"*L"x^4"*"$(λ_coefs[4])"*L"x^3"*"+$(λ_coefs[3])"*L"x^2"*"$(λ_coefs[2])"*L"x"*"$(λ_coefs[1])")
push!(AllPlots,plt_λLN);    push!(plt_img_names,"lambda_ATr_poly")

#############################################################################################
# Distributions of A̅ by binning joint probability distribution P(t̅ᵣ,A̅)
t̅ᵣ_bins, A̅_WeiFits, A̅_WeiModes, plt_tcbins, plt_BinDist = distros_from_jpd(allt̅ᶜᵣ, allA̅, Weibull,L"\overline{t}_d")
plot!(plt_tcbins, xlab = L"\overline{t}_d", ylab = L"\overline{A}", title=L"P(\overline{t}_d,\overline{A})",xlim=(-t̅ₑ,t̅ₑ), ylim=(minimum(allA̅),maximum(allA̅)))
plot!(plt_BinDist, xlab=L"\overline{A}", ylab=L"P(\overline{A})", palette=:darkrainbow)
push!(AllPlots,plt_BinDist);    push!(plt_img_names,"Atr_BinDists")

# Interpolation of shape and scale parameters of resulting Weibull distributions per bin
Nᴮ = length(t̅ᵣ_bins)
## B-spline interpolation for κ
κ_itp = interpolate(t̅ᵣ_bins, A̅_WeiFits[:,1], BSplineOrder(4))    
# t̅ᵣ_range = range(t̅ᵣˡᵇ+wᴮ/2,t̅ᵣᵘᵇ-wᴮ/2,1000)
t̅ᵣ_range = range(t̅ᵣ_bins[1],t̅ᵣ_bins[end],1000)
κⁱ = κ_itp.(t̅ᵣ_range)  

plt_kw = plot(title="A̅ Distributions: Weibull shape parameter", xlab=L"\overline{t}_d", ylab=L"\kappa", palette=[cb[11];cb[8];cb[4]])
plot!(plt_kw, t̅ᵣ_bins, A̅_WeiFits[:,1], seriestype=:scatter, lab=L"\kappa"*"(bins)")
plot!(plt_kw, t̅ᵣ_range, κⁱ, lab=L"\kappa (\overline{t}_d)"*" B-Spline interpolation", lw=2)

## Interpolation with Gaussian function
# sₜ =  t̅ᵣ_bins[Int(floor(Nᴮ/2))+1] # x-axis shift 
sₜ =  t̅ᵣ_range[findmax(κⁱ)[2]] # x-axis shift 
sₖ = median(A̅_WeiFits[:,1])     # y-axis shift
# K =  A̅_WeiFits[Int(floor(Nᴮ/2))+1,1] - sₖ   # Peak value of function
K =  maximum(κⁱ) - sₖ   # Peak value of function
# Dispersion
p1 = findall(diff(sign.(κⁱ.-sₖ.-K*exp(-1/4))) .== 2)
p2 = findall(diff(sign.(κⁱ.-sₖ.-K*exp(-1/4))) .== -2)
Tₖ = abs(t̅ᵣ_range[p2[1]] - t̅ᵣ_range[p1[1]])  

gᵏ = cgaussian(t̅ᵣ_range, K, sₜ, Tₖ) .+ sₖ
plot!(plt_kw, t̅ᵣ_range, gᵏ, lw=2, lab=L"\kappa (\overline{t}_d)=K \exp{(\frac{\overline{t}_d-s_t}{T_{\kappa}}})^2+s_{\kappa}")
push!(AllPlots,plt_kw);    push!(plt_img_names,"k_Weib_interp")

## B-spline interpolation for λ
λ_itp = interpolate(t̅ᵣ_bins, A̅_WeiFits[:,2], BSplineOrder(4))    
λⁱ = λ_itp.(t̅ᵣ_range)  

plt_λw = plot(title="A̅ Distributions: Weibull scale parameter", xlab=L"\overline{t}_d", ylab=L"\lambda", palette=[cb[11];cb[8];cb[4]])
plot!(plt_λw, t̅ᵣ_bins, A̅_WeiFits[:,2], seriestype=:scatter, lab=L"\lambda"*"(bins)")
plot!(plt_λw, t̅ᵣ_range, λⁱ, lab=L"\lambda (\overline{t}_d)"*" B-Spline interpolation", lw=2)

## Interpolation with Gaussian function
sₜ =  t̅ᵣ_range[findmax(λⁱ)[2]] # x-axis shift 
sₗ = median(A̅_WeiFits[:,2]) # y-axis shift 
# Λ =  A̅_WeiFits[Int(floor(Nᴮ/2))+1,2] - sₗ # Peak value of function
Λ =  maximum(λⁱ) - sₗ # Peak value of function
p1 = findall(diff(sign.(λⁱ.-sₗ.-Λ*exp(-1/4))) .== 2)
p2 = findall(diff(sign.(λⁱ.-sₗ.-Λ*exp(-1/4))) .== -2)
# Dispersion
Tₗ = abs(t̅ᵣ_range[p2[1]] - t̅ᵣ_range[p1[1]])

gˡ = cgaussian(t̅ᵣ_range, Λ, sₜ, Tₗ) .+ sₗ
plot!(plt_λw, t̅ᵣ_range, gˡ, lw=2, lab=L"\lambda (\overline{t}_d)=\Lambda \exp{(\frac{\overline{t}_d-s_t}{T_{\lambda}}})^2+s_{\lambda}")
push!(AllPlots,plt_λw);    push!(plt_img_names,"λ_Weib_interp")

#############################################################################################
# Distributions of T̅ᵣ by binning joint probability distribution P(A̅,T̅ᵣ)
A̅_bins, T̅ᵣ_LogNormFits, T̅ᵣ_LogNormModes, plt_Abins, plt_BinDist2 = distros_from_jpd(allA̅, allT̅ᵣ, LogNormal,L"\overline{A}")
plot!(plt_Abins, ylab = L"\overline{T}_d", xlab = L"\overline{A}", title=L"P(\overline{A},\overline{T}_d)", ylim=(minimum(allT̅ᵣ),1), xlim=(minimum(allA̅),1))
plot!(plt_BinDist2, xlab=L"\overline{T}_d", ylab=L"P(\overline{T}_d)", palette=:darkrainbow)
push!(AllPlots,plt_BinDist2);    push!(plt_img_names,"Tr_LogNorm_BinDists")

# Interpolation of shape and scale parameters of resulting Weibull distributions per bin
## B-spline interpolation for κ
μ_itp = interpolate(A̅_bins, T̅ᵣ_LogNormFits[:,1], BSplineOrder(4))    
A̅_range = range(A̅_bins[1],A̅_bins[end],1000)
μⁱ = μ_itp.(A̅_range)  

plt_μLN = plot(title="T̅ᵣ Distributions: LogNormal shape parameter", xlab=L"\overline{A}", ylab=L"\mu", palette=[cb[11];cb[8];cb[4]])
plot!(plt_μLN, A̅_bins, T̅ᵣ_LogNormFits[:,1], seriestype=:scatter, lab=L"\mu"*"(bins)")
plot!(plt_μLN, A̅_range, μⁱ, lab=L"\mu (\overline{A})"*" B-Spline interpolation", lw=2)

## Polynomial interpolation
μ_poly = pfit(A̅_bins, T̅ᵣ_LogNormFits[:,1],3); μ_coefs = round.(μ_poly[:]*1000)./1000
μ_eval = μ_poly.(A̅_range)

plot!(plt_μLN, A̅_range, μ_eval, lw=2, lab=L"\mu (\overline{A})"*"=$(μ_coefs[4])"*L"~x^3"*"+$(μ_coefs[3])"*L"~x^2"*"$(μ_coefs[2])"*L"~x"*"$(μ_coefs[1])")
push!(AllPlots,plt_μLN);    push!(plt_img_names,"mu_Tr_poly")

σ_itp = interpolate(A̅_bins, T̅ᵣ_LogNormFits[:,2], BSplineOrder(4))    
σⁱ = σ_itp.(A̅_range)  

plt_σLN = plot(title="T̅ᵣ Distributions: LogNormal scale parameter", xlab=L"\overline{A}", ylab=L"\sigma", palette=[cb[11];cb[8];cb[4]])
plot!(plt_σLN, A̅_bins, T̅ᵣ_LogNormFits[:,2], seriestype=:scatter, lab=L"\sigma"*"(bins)")
plot!(plt_σLN, A̅_range, σⁱ, lab=L"\sigma (\overline{A})"*" B-Spline interpolation", lw=2)

## Polynomial interpolation
σ_poly = pfit(A̅_bins, T̅ᵣ_LogNormFits[:,2],3); σ_coefs = round.(σ_poly[:]*1000)./1000
σ_eval = σ_poly.(A̅_range)

plot!(plt_σLN, A̅_range, σ_eval, lw=2, lab=L"\sigma (\overline{A})"*"=$(σ_coefs[4])"*L"~x^3"*"+$(σ_coefs[3])"*L"~x^2"*"$(σ_coefs[2])"*L"~x"*"+$(σ_coefs[1])")
push!(AllPlots,plt_σLN);    push!(plt_img_names,"sigma_Tr_poly")

#############################################################################################
# Distributions of T̅ᵣ by binning joint probability distribution P(t̅ᵣ,T̅ᵣ)
t̅ᵣ_bins, T̅ᵣ_LogNormFits, T̅ᵣ_LogNormModes, plt_bins_trTr, plt_BinDist3 = distros_from_jpd(allt̅ᶜᵣ, allT̅ᵣ, LogNormal,L"\overline{t}_d")

plot!(plt_BinDist3, xlab=L"\overline{T}_d", ylab=L"P(\overline{T}_d)", palette=:darkrainbow)
push!(AllPlots,plt_BinDist3);    push!(plt_img_names,"Trtr_LogNorm_BinDists")

## B-spline interpolation for μ
μ_itp = interpolate(t̅ᵣ_bins, T̅ᵣ_LogNormFits[:,1], BSplineOrder(4))    
t̅ᵣ_range = range(t̅ᵣ_bins[1],t̅ᵣ_bins[end],1000)
μⁱ = μ_itp.(t̅ᵣ_range)  

plt_μLN = plot(title="T̅ᵣ Distributions: LogNormal shape parameter", xlab=L"\overline{t}_d", ylab=L"\mu", palette=[cb[11];cb[8];cb[4]])
plot!(plt_μLN, t̅ᵣ_bins, T̅ᵣ_LogNormFits[:,1], seriestype=:scatter, lab=L"\mu"*"(bins)")
plot!(plt_μLN, t̅ᵣ_range, μⁱ, lab=L"\mu (\overline{t}_d)"*" B-Spline interpolation", lw=2)

## Polynomial interpolation
μ_poly_t̅ᵣT̅ᵣ = pfit(t̅ᵣ_bins, T̅ᵣ_LogNormFits[:,1],8); μ_coefs = round.(μ_poly_t̅ᵣT̅ᵣ[:]*1000)./1000
μ_eval = μ_poly_t̅ᵣT̅ᵣ.(t̅ᵣ_range)

plot!(plt_μLN, t̅ᵣ_range, μ_eval, lw=2)#, lab=L"\mu (\overline{t}_d)"*"=$(μ_coefs[5])"*L"x^4"*"$(μ_coefs[4])"*L"x^3"*"+$(μ_coefs[3])"*L"x^2"*"$(μ_coefs[2])"*L"x"*"$(μ_coefs[1])")
push!(AllPlots,plt_μLN);    push!(plt_img_names,"mu_Trtr_poly")

## B-spline interpolation for σ
σ_itp = interpolate(t̅ᵣ_bins, T̅ᵣ_LogNormFits[:,2], BSplineOrder(4))    
σⁱ = σ_itp.(t̅ᵣ_range)  

plt_σLN = plot(title="T̅ᵣ Distributions: LogNormal scale parameter", xlab=L"\overline{t}_d", ylab=L"\sigma", palette=[cb[11];cb[8];cb[4]])
plot!(plt_σLN, t̅ᵣ_bins, T̅ᵣ_LogNormFits[:,2], seriestype=:scatter, lab=L"\sigma"*"(bins)")
plot!(plt_σLN, t̅ᵣ_range, σⁱ, lab=L"\sigma (\overline{t}_d)"*" B-Spline interpolation", lw=2)

## Polynomial interpolation
σ_poly_t̅ᵣT̅ᵣ = pfit(t̅ᵣ_bins, T̅ᵣ_LogNormFits[:,2],9); σ_coefs = round.(σ_poly_t̅ᵣT̅ᵣ[:]*1000)./1000
σ_eval = σ_poly_t̅ᵣT̅ᵣ.(t̅ᵣ_range)

plot!(plt_σLN, t̅ᵣ_range, σ_eval, lw=2)#, lab=L"\sigma (\overline{t}_d)"*"=$(σ_coefs[5])"*L"x^4"*"+$(σ_coefs[4])"*L"x^3"*"$(σ_coefs[3])"*L"x^2"*"$(σ_coefs[2])"*L"x"*"+$(σ_coefs[1])")
push!(AllPlots,plt_σLN);    push!(plt_img_names,"sigma_Trtr_poly")

#############################################################################################
# Sampling of A,T and tᶜ parameters for NoGEEs number of Gaussian Elementary Envelopes 
T̅ᵣ, A̅ = [Ø[:] for _ = 1:2]
u = zeros(Float64,Nₜ)
pltCONV = plot()

# t̅ᶜᵣ = -median(t̅ᶜᵣ_Fit)
t̅ᶜᵣ = [rand(Δtᶜᵣ_Fit)]
wΔ = rand(Δtᶜᵣ_Fit,NoGEEs-1)
for i ∈ 1:Int((NoGEEs-1)/2)
    Σₗ = t̅ᶜᵣ[1]-sum(wΔ[1:i])
    Σᵣ = t̅ᶜᵣ[1] + sum(wΔ[Int((NoGEEs-1)/2)+1:Int((NoGEEs-1)/2)+i])
    push!(t̅ᶜᵣ,Σₗ)
    push!(t̅ᶜᵣ,Σᵣ)
end
t̅ᶜᵣ = sort(t̅ᶜᵣ)
t̅ᶜ = t̅ᶜᵣ.*(T₀ₚ./Tₚⱼ)
tᶜ = t̅ᶜ.*Tₚⱼ
t̅ᶜₑ = tᶜ/(ΔTₑᵥᴴ * ϵₛₛ*T₀ₛ)

#############################################################################################
# ## Determination of T
# T̅ᵣ_range = range(0, maxT̅ᵣ, length=2^9)
# plt_LogNorms = plot(xlab=L"\overline{T}_d" ,ylab="Density", palette=:darkrainbow)
# for i ∈ 1:NoGEEs
#     # A: Sampling from dataset by binning using interpolated LogNormal parameters
#     μᴸᴺ = μ_poly_t̅ᵣT̅ᵣ(t̅ᶜᵣ[i])
#     σᴸᴺ = σ_poly_t̅ᵣT̅ᵣ(t̅ᶜᵣ[i])
#     fitted_dist = LogNormal(μᴸᴺ,σᴸᴺ)

#     pdf_T̅ᵣ = pdf.(fitted_dist, T̅ᵣ_range)

#     T̅ᵣⁱ = rand(fitted_dist)
#     while T̅ᵣⁱ < Tᶜᵘᵗ/T₀ₚ
#         T̅ᵣⁱ = rand(fitted_dist)
#     end
#     push!(T̅ᵣ, T̅ᵣⁱ)
#     plot!(plt_LogNorms, T̅ᵣ_range,pdf_T̅ᵣ, lw=2, lab=L"P(\overline{T}_d)"*"for GEE $i")
# end
# push!(AllPlots,plt_LogNorms);    push!(plt_img_names,"Tr_LogNorms_sampled")

#############################################################################################
# Conditional sampling for A̅
wᴮ = 2*std(Δtᶜᵣ_Fit)  # For A1, B
# wᴮ = 2*std(Δtᶜᵣ_Fit)*T₀ₚ/(ΔTₑᵥᴴ * ϵₛₛ*T₀ₛ)  # For A2
A̅_range = range(0, maxA̅, length=2^9)
# plt_weibs = plot3d(xlab=L"\overbar{t}^c", ylab=L"\overbar{A}" ,zlab="Density")
plt_weibs = plot(xlab=L"\overline{A}" ,ylab="Density", palette=:darkrainbow)

for i ∈ 1:NoGEEs
    # A1: Sampling from copula 
    # _, A̅_samp, _ = conditional_sampling(JPD_tᶜA,10000,t̅ᶜ[i],wᴮ)
    # pdf_A̅, fitted_dist, _ = dist_fit(A̅_samp, Weibull, 2^9)

    # A2: Sampling from copula 
    # _, A̅_samp, _ = conditional_sampling(JPD_tᶜₑA,10000,t̅ᶜₑ[i],wᴮ)
    # pdf_A̅, fitted_dist, _ = dist_fit(A̅_samp, Weibull, 2^9)

    # B1: Sampling from dataset by binning using interpolated Weibull parameters
    κ = K*exp(-((t̅ᶜᵣ[i]-sₜ)/Tₖ)^2) + sₖ
    λ = Λ*exp(-((t̅ᶜᵣ[i]-sₜ)/Tₗ)^2) + sₗ
    fitted_dist = Weibull(κ,λ)
    pdf_A̅ = pdf.(fitted_dist, A̅_range)

    # # B2: Sampling from dataset by binning using interpolated Weibull parameters
    # κ = k_poly(T̅ᵣ[i])
    # λ = λ_poly(T̅ᵣ[i])
    # fitted_dist = Weibull(κ,λ)
    # pdf_A̅ = pdf.(fitted_dist, A̅_range)

    push!(A̅, rand(fitted_dist))
    # plot3d!(plt_weibs, ones(2^9)*t̅ᶜᵣ[i], A̅_range,pdf_A̅, lw=2, lab=L"P(A̅)"*"for GEE $i")
    plot!(plt_weibs, A̅_range,pdf_A̅, lw=2, lab=L"P(\overline{A})"*"for GEE $i")
end
push!(AllPlots,plt_weibs);    push!(plt_img_names,"A_Weibs_sampled")

# Dimensionalise based on given Hₛ, Tₚ
A = A̅.*Hₛⱼ

#############################################################################################
## Determination of T
# Conditional sampling for T̅ᵣ
# wᴮ = std(A̅_Fit)/2     # For A1
wᴮ = 2*std(Δtᶜᵣ_Fit)*T₀ₚ/(ΔTₑᵥᴴ * ϵₛₛ*T₀ₛ)  # For A2
T̅ᵣ_range = range(0, maxT̅ᵣ, length=2^9)
plt_LogNorms = plot(xlab=L"\overline{T}_d" ,ylab="Density", palette=:darkrainbow)

for i ∈ 1:NoGEEs
    # B1: Sampling from (A̅,T̅ᵣ) copula
    # _, T̅ᵣ_samp,_ = conditional_sampling(JPD_ATᵣ,1000,A̅[i],wᴮ)
    # pdf_T̅ᵣ, fitted_dist, _ = dist_fit(T̅ᵣ_samp, LogNormal, 2^9)
    
    # B2: Sampling from (t̅ᵣ,T̅ᵣ) copula
    # _, T̅ᵣ_samp,_ = conditional_sampling(JPD_tᶜₑTᵣ,10000,t̅ᶜₑ[i],wᴮ)
    # pdf_T̅ᵣ, fitted_dist, _ = dist_fit(T̅ᵣ_samp, LogNormal, 2^9)

    # A: Sampling from dataset by binning using interpolated LogNormal parameters
    μᴸᴺ = μ_poly(A̅[i])
    σᴸᴺ = σ_poly(A̅[i])
    fitted_dist = LogNormal(μᴸᴺ,σᴸᴺ)

    pdf_T̅ᵣ = pdf.(fitted_dist, T̅ᵣ_range)

    T̅ᵣⁱ = rand(fitted_dist)
    while T̅ᵣⁱ < Tᶜᵘᵗ/T₀ₚ
        T̅ᵣⁱ = rand(fitted_dist)
    end
    push!(T̅ᵣ, T̅ᵣⁱ)
    # plot3d!(plt_weibs2, ones(2^9)*A̅[i], T̅ᵣ_range,pdf_T̅ᵣ, lw=2, lab=L"P(T̅ᵣ)"*"for GEE $i")
    plot!(plt_LogNorms, T̅ᵣ_range,pdf_T̅ᵣ, lw=2, lab=L"P(\overline{T}_d)"*"for GEE $i")
end
push!(AllPlots,plt_LogNorms);    push!(plt_img_names,"Tr_LogNorms_sampled")

# Dimensionalise based on given Tₚ
T̅ = T̅ᵣ.*(T₀ₚ./Tₚⱼ)
T = T̅ .* Tₚⱼ
G = gauss_fun(t, A, tᶜ, T)

#############################################################################################
# Plot sampled values on top of corresponding joint pdfs
# push!(AllPlots,hplt_TᵣA);    push!(plt_img_names,"hmap_P(Td,A)")
# push!(AllPlots,hplt_tᶜₑA);    push!(plt_img_names,"hmap_P(tev,A)")
# push!(AllPlots,plt_tcbins);    push!(plt_img_names,"tc_bins")
# push!(AllPlots,plt_Abins);    push!(plt_img_names,"A_Bins")

plot!(hplt_TᵣA, [T̅ᵣ], [A̅], seriestype=:scatter, lab="Sampled")
push!(AllPlots,hplt_TᵣA);    push!(plt_img_names,"TrA_samp_on_hist")
plot!(hplt_tᶜₑA, [t̅ᶜₑ], [A̅], xlim=(-1,1), ylim=(minimum(allA̅),maximum(allA̅)), seriestype=:scatter, lab="Sampled")
push!(AllPlots,hplt_tᶜₑA);    push!(plt_img_names,"tceA_samp_on_hist")
plot!(plt_tcbins, t̅ᶜᵣ, A̅, seriestype=:scatter, lab="Sampled", color=:green)
push!(AllPlots,plt_tcbins);    push!(plt_img_names,"tc_bins_wsamp")
plot!(plt_Abins, A̅, T̅ᵣ, seriestype=:scatter, lab="Sampled")
push!(AllPlots,plt_Abins);    push!(plt_img_names,"A_Bins_wsamp")

#############################################################################################
# Reconstruction
Nₛ₂ = Int(nextpow(2,Nₜ)/2)
TOT = zeros(ComplexF64, Nₜ)  # Total surface elevation
SP_TOT = zeros(Float64, Nₛ₂) # Total spectrum of propagated Gaussian EWGs

plt_gn = plot(palette=:darkrainbow)
pltEWG = plot(palette=:darkrainbow)
pltSPEC = plot(palette=:darkrainbow)
plt_demo = plot(palette=[cb[8];cb[11];cb[4]])

prop = 0    # Reconstruction method
# Random sampling of β
# β̃ = fitβ[2]*tᶜ    # Directly from linear fit
β̃ = rand(β̇_Fit)*t̅ᶜᵣ # Through β̇
plt_tᶜβ = plot(tᶜ, fitβ[2]*tᶜ, xlab = L"t~[s]", ylab = L"\beta~[rad]", lab="Linear fit")
plot!(plt_tᶜβ, tᶜ, β̃, seriestype=:scatter, lab="Sampled values")
push!(AllPlots,plt_tᶜβ);    push!(plt_img_names,"tcbeta_linfit")
# Random sampling of omega
Ω = rand(Ω_Fit)
# Ωₙ = rand(Ω_Fit,NoGEEs)
# Additional parameters
ωₚ = 2π/Tₚⱼ
ωᵧ = 4*acos(exp(-1/4))./T
ω₁ = Ω_dist[1]*ωₚ
ω₂ = 0.0

for n ∈ 1:NoGEEs
    # Propagated EWG
    global ω₂ = ω₁ - ωᵧ[n]
    gₙ, ηᶜ, ηₙ, FR, MAG = recon(prop, A̅[n], t̅ᶜ[n], T̅[n], β̃[n], Hₛⱼ, Tₚⱼ, Ω, t)
    # gₙ, ηᶜ, ηₙ, FR, MAG = recon(prop, A̅[n], t̅ᶜ[n], T̅[n], β̃[n], Hₛⱼ, Tₚⱼ, Ωₙ[n], t)

    TOT[:] = TOT[:] .+ real.(ηₙ)
    SP_TOT[:] = SP_TOT[:] .+ MAG

    plot!(plt_gn, t, gₙ, xlab = "t [s]", ylab = "A [m]", lab = "g$(n)", lw=2, legend=:outerbottom, legendcolumns=min(NoGEEs,8))
    plot!(pltEWG, t, real(ηₙ), xlab = "t [s]", ylab = "η [m]", lab = "η$(n)", lw=2, legend=:outerbottom, legendcolumns=min(NoGEEs,8))
    plot!(pltSPEC, FR, MAG, xlab = L"f [Hz]", ylab = L"A [m]", lab = "g$(n)", lw=2, legend=:outerbottom, legendcolumns=min(NoGEEs,8))
    plot!(pltSPEC, xlim=(0,fᶜᵘᵗ))

    if n == 1
        plot!(plt_demo, t, gₙ, xlab = "t [s]", ylab = "[m]", lab = L"g_1", lw=2)
        plot!(plt_demo, t, real(ηᶜ), lab="Temporal term", lw=2)
        plot!(plt_demo, t, real(ηₙ), lab="Propagated elementary envelope", lw=2)
        # plot!(plt_demo, xlim=(-4*T̅ₒ[n]*Tₚ,4*T̅ₒ[n]*Tₚ))
    end                         
end

push!(AllPlots,plt_gn);    push!(plt_img_names,"EGEs")
push!(AllPlots,pltEWG);    push!(plt_img_names,"GFWGs")
push!(AllPlots,pltSPEC);    push!(plt_img_names,"Gspectra")

ηᴳ = real(TOT) 
plt_etaG = plot(t,ηᴳ,xlab="t [s]", ylab="η [m]", lw=2, lab="HyperGroup") #, xlim=(t[Int(Nₜ/2-2^10)], t[Int(Nₜ/2+2^10)]))
push!(AllPlots,plt_etaG);    push!(plt_img_names,"HyperGroup")

uᴳ, u̇ᴳ = derivatives(ηᴳ,t)
plt_udotG = plot(t,u̇ᴳ,xlab="t [s]", ylab="u̇ [m/s²]", lw=2, lab="acceleration")

fᴳ, Sₐᴳ, θᴳ, df,_ = one_side_asp(ηᴳ ,t)
Sᴳ = 0.5/df .* Sₐᴳ.^2
Eᴳ = sum(Sᴳ*df)     # Spectral energy
Hₛᴳ = 4*std(ηᴳ)
Tₚᴳ = 1/fᴳ[findmax(Sᴳ)[2]]
H̅ₛᴳ = Hₛᴳ/Hₛⱼ  # Normalised significant wave height
T̅ₚᴳ =  Tₚᴳ/Tₚⱼ   # Normalised peak wave period
ϵₛₛᴳ = Hₛᴳ/Tₚᴳ^2 * 2π/g  # Effective steepness

# dfⱼ = (fᶜᵘᵗ-1/Tₑ)/(Nₛ₂-1)
# fⱼ,Sⱼ = spectrum(Hₛⱼ,Tₚⱼ,3.3,Tᶜᵘᵗ,Tₑ,Nₛ₂)
# Eⱼ = sum(Sⱼ*dfⱼ) 
# Eⱼ = Hₛⱼ^2 / 16
fⱼ,Sⱼ = spectrum(Hₛᴳ,Tₚⱼ,3.3,Tᶜᵘᵗ,Tₑ,Nₛ₂)
Eⱼ = Hₛᴳ^2 / 16

plt_HGS = plot(tile="Spectral comparison - Tₚ=$Tₚⱼ", xlab = L"f [Hz]", ylab = L"S(f)~[m^2 s]", palette=[cb[11];cb[4];cb[8]])
plot!(plt_HGS, fᴳ, Sᴳ, lab = L"S_G:~H_s="*"$(round(4*std(ηᴳ)*1000)/1000) , E=$(round(Eᴳ*1000)/1000)", xlim=(0, fᶜᵘᵗ), lw=2)
# plot!(plt_HGS, fⱼ, Sⱼ, lab = L"S_j:~H_s="*"$Hₛⱼ , E=$(round(Eⱼ*1000)/1000)", xlim=(0, fᶜᵘᵗ), lw=2, ls=:dash)
plot!(plt_HGS, fⱼ, Sⱼ, lab = L"S_j:~H_s="*"$(round(4*std(ηᴳ)*1000)/1000) , E=$(round(Eⱼ*1000)/1000)", xlim=(0, fᶜᵘᵗ), lw=2, ls=:dash)
push!(AllPlots,plt_HGS);    push!(plt_img_names,"HG_spectrum")

# Gaussian approximation against envelope and scatter plot of peaks
plt_envs = plot(xlab = L"t~[s]", ylab = L"[m]", palette=[cb[4];cb[8];cb[11]], legend=:outertop, legendcolumns=4)

if gen_method == 1
    plot!(plt_envs,t, [u G], lab = [L"u(t)" L"G(t)"], lw=[2 2])
else
    # plot!(plt_envs,t, [u G], lab = [L"u(t)" L"G(t)"], lw=[2 2])
    plot!(t, G, lab = L"G(t)", lw=2)
end
plot!(plt_envs,tᶜ, A, seriestype=:scatter, ms=2, mc=:red, lab = "Peaks")
push!(AllPlots,plt_envs);    push!(plt_img_names,"Envelope")

#############################################################################################
# Write Outputs
if wout
    wpath = joinpath(OFASTpath,"ExtElev",case_id)
    # List all files in the directory
    files_ExtElev = readdir(wpath)
    # Filter files that start with "DWE_" and count them
    HG_fcount = count(file -> startswith(file, "HG_"), files_ExtElev)
    # Write the HyperGroup file
    fid = joinpath(wpath,"HG_$(HG_fcount+1).Elev")
    open(fid, "w")
    writedlm(fid, [t.+t[end] ηᴳ], '\t')

    fid = joinpath(stats_path,"HyperGroup","HG_$(HG_fcount+1)_parameters")
    head = ["t̅ᶜᵣ [-]" "T̅ᵣ [-]" "A̅ [-]" "β̃ [rad]" "" ""]
    head2 = ["Ω" "T̅ₚᴳ" "H̅ₛᴳ" "Eᴳ" "ϵₛₛᴳ" "ΔTₑᵥᴳ"]
    row = round.([t̅ᶜᵣ T̅ᵣ A̅ β̃]*1e3)./1e3
    ∅ₛ = ["" for _ in 1:NoGEEs]
    empty_row = ["" for _ in 1:6]
    row = [row ∅ₛ ∅ₛ]
    row2 = round.([Ω T̅ₚᴳ H̅ₛᴳ Eᴳ ϵₛₛᴳ ΔTₑᵥᴴ]*1e3)./1e3
    open(fid, "w")
    writedlm(fid, [head; row; reshape(empty_row,1,:); head2; row2], '\t')
end

#############################################################################################
# Plots
for i ∈ 1:length(AllPlots)
    display(AllPlots[i])
    savefig(AllPlots[i],joinpath(stats_path,"HyperGroup","$(plt_img_names[i]).png"))
    savefig(AllPlots[i],joinpath(stats_path,"HyperGroup","$(plt_img_names[i]).svg"))
end

end

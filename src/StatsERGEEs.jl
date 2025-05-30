# Load necessary modules
using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles
using FFTW, DSP, Statistics, BSplineKit, StatsBase, Distributions, Copulas
using CurveFit: linear_fit
using StatsPlots
import ColorSchemes.darkrainbow

gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))
cb = darkrainbow

# Include necessary scripts for functions
include("func/text_process.jl"); include("func/signal_processing.jl")
include("func/peak_detect.jl");  include("func/runge_kutta.jl")
include("func/fugees.jl");  include("func/directories.jl")
include("func/wave_theory.jl");  include("func/jonswap.jl")

#############################################################################################
const ρ, g, d::Float64 = 1025.0, 9.81, 100.0  # density [kg/m³], gravity [m/s²], depth [m]
Tₚⱼ = 10.0
Hₛⱼ = 8.0
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

MaxNoEv = maximum(NoEv)
allG = zeros(Float64,Nₜ,TotNoEv)
all_sp2, all_phi = zeros(Float64,Nₛ₂,TotNoEv), zeros(Float64,Nₛ₂,TotNoEv)
fsp2 = range(0.0,4.0,Nₛ₂)
allΩ, allωₚ, allβ, alltᶜ, allAₒ, allTₒ, allT₂₋, allΔTₑᵥ, allΔtᶜ, maxAₒ  = [], [], [], [], [], [], [], [], [], []
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
        push!(allΩ, Ω)
        push!(allωₚ, 2π/Tₚ)
        push!(maxAₒ, maximum(A̅ₒ))

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

t = t̅*T̃ₚᴿ
# η = real(Hₛⱼ*Gmean.*exp.(-1im * mean(allΩ)*2π/Tₚⱼ * t̅*Tₚⱼ))

plt_eta2 = plot(t, η₂₋, xlab = "t [s]", ylab = L"\eta_{2nd} ~[m]", lw=2)

# STATISTICAL ANALYSIS
allΩ = Float64.(allΩ);  allβ = Float64.(allβ) 
allAₒ = Float64.(allAₒ);    allTₒ = Float64.(allTₒ);    alltᶜ = Float64.(alltᶜ)
allΔtᶜ = Float64.(allΔtᶜ);  allT₂₋ = Float64.(allT₂₋);  allΔTₑᵥ = Float64.(allΔTₑᵥ) 

## Scatter plots
pltΩ = plot([allΩ], seriestype=:scatter, legend=:false ,xlab="Event", ylab=L"\Omega=\frac{\tilde{\omega}}{\omega_p}~[-]")
pltT2 = plot([allT₂₋], seriestype=:scatter, legend=:false ,xlab="Event", ylab=L"T_{2^-}~[s]")
pltΔT = plot([allΔTₑᵥ], seriestype=:scatter, legend=:false ,xlab="Event", ylab=L"\Delta T_{ev}~[s]")
pltFair = plot([FairTen], seriestype=:scatter, legend=:false ,xlab="Event", ylab="Fairlead Tension [kN]")

## Histograms
### Ω
μΩ = mean(allΩ);    σΩ = std(allΩ); 
Nrang = range(μΩ-4*σΩ, μΩ+4*σΩ, 2^9);  Ncurve = 1/(σΩ*sqrt(2π)) * exp.(-((Nrang.-μΩ)/(sqrt(2)*σΩ)).^2)
hplt_Ω = histogram(allΩ, normalize=:pdf, lab=:false)
plot!(hplt_Ω, title="Probability density function: P(Ω)", xlab=L"Ω", ylab=L"P(Ω)")
plot!(hplt_Ω, Nrang, Ncurve, lw=2, lab="Normal Fit: μ=$(round(μΩ*1e3)/1e3), σ=$(round(σΩ*1e3)/1e3)")

### Δtᶜ
WeiΔtᶜ = fit(Weibull, allΔtᶜ);  kʷ = WeiΔtᶜ.α;    λʷ = WeiΔtᶜ.θ
Wrang = range(0, maximum(allΔtᶜ), 2^9);  Wcurve = kʷ/λʷ * (Wrang/λʷ).^(kʷ-1) .* exp.(-(Wrang/λʷ).^kʷ)
hplt_Dtc = histogram(allΔtᶜ, normalize=:pdf, lab=:false)
plot!(hplt_Dtc, title="Probability density function: P(Δtᶜ)", xlab=L"\Delta t^c", ylab=L"P(\Delta t^c)")
plot!(hplt_Dtc, Wrang, Wcurve, lw=2, lab="Weibull Fit: κ=$(round(kʷ*1e3)/1e3), λ=$(round(λʷ*1e3)/1e3)")

### Tₒ
WeiTₒ = fit(Weibull, allTₒ);  kʷ = WeiTₒ.α;    λʷ = WeiTₒ.θ
Wrang = range(0, maximum(allTₒ), 2^9);  Wcurve = kʷ/λʷ * (Wrang/λʷ).^(kʷ-1) .* exp.(-(Wrang/λʷ).^kʷ)
hplt_To = histogram(allTₒ, normalize=:pdf, lab=:false)
plot!(hplt_To, title="Probability density function: P(Tₒ)", xlab=L"T_o", ylab=L"P(T_o)")
plot!(hplt_To, Wrang, Wcurve, lw=2, lab="Weibull Fit: κ=$(round(kʷ*1e3)/1e3), λ=$(round(λʷ*1e3)/1e3)")

### Aₒ
WeiAₒ = fit(Weibull, allAₒ);  kʷ = WeiAₒ.α;    λʷ = WeiAₒ.θ
Wrang = range(0, maximum(allAₒ), 2^9);  Wcurve = kʷ/λʷ * (Wrang/λʷ).^(kʷ-1) .* exp.(-(Wrang/λʷ).^kʷ)
hplt_Ao = histogram(allAₒ, normalize=:pdf, lab=:false)
plot!(hplt_Ao, title="Probability density function: P(Aₒ)", xlab=L"A_o", ylab=L"P(A_o)")
plot!(hplt_Ao, Wrang, Wcurve, lw=2, lab="Weibull Fit: κ=$(round(kʷ*1e3)/1e3), λ=$(round(λʷ*1e3)/1e3)")

# Copulas
AₒCDF = cdf.(WeiAₒ,allAₒ)
TₒCDF = cdf.(WeiTₒ,allTₒ)
# AT_cop = fit(FrankCopula,[TₒCDF AₒCDF])
# # Create a grid for plotting
# grid_size = length(AT_cop)
# A_grid = collect(range(0, stop=1, length=grid_size))
# T_grid = collect(range(0, stop=1, length=grid_size))
# # Evaluate the copula density on the grid
# dens = zeros(grid_size, grid_size)
# dens = [pdf(AT_cop, [Ti, Ai]) for Ti in T_grid, Ai in A_grid]
# for (i, Ti) in enumerate(T_grid)
#     for (j, Ai) in enumerate(A_grid)
#         dens[i, j] = pdf(AT_cop, [Ti; Ai])
#     end
# end
# # Plot the copula density
# contour(T_grid, A_grid, dens, fill=true, c=:viridis, title="Gaussian Copula Density", xlabel="To ", ylabel="Ao")

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

### tᶜ-Tₒ
hplt_tcT = histogram2d(alltᶜ, allTₒ, normalize=:pdf, bins=(tcbins,Tobins), show_empty_bins=true, color=:plasma, xlab = L"\overline{t}", ylab = L"\overline{T}")

Aint = zeros(Float64, size(probabs)[2])

for i ∈ 1:size(probabs)[2]
    j = findmax(probabs[:,i])[2]
    Aint[i] = eA[j]
end

A_itp = interpolate(etc[1:end-1], Aint, BSplineOrder(4))
Aᵢ = A_itp.(t̅)

Aₒᵢ = zeros(Float64,Nₜ)
Aₒᵢ[Int(Nₜ/2-2^9):Int(Nₜ/2+2^9)] = Aᵢ[Int(Nₜ/2-2^9):Int(Nₜ/2+2^9)]

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
    u = Aₒᵢ*H̃ₛᴿ
    wdtc = rand(WeiΔtᶜ,6)
    tᶜₒ = [-sum(wdtc[1:3]); -sum(wdtc[1:2]); -wdtc[1]; 0; wdtc[4]; sum(wdtc[4:5]); sum(wdtc[4:6])].*T̃ₚᴿ
    i⁺ = Int.(round.(tᶜₒ./(t[2]-t[1]) .+ Int(Nₜ/2)))
    Aₒ = u[i⁺]
    Tₒ = rand(WeiTₒ,7).*T̃ₚᴿ
    G = gauss_fun(t, Aₒ, tᶜₒ, Tₒ)
end

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

NDΩ = Normal(μΩ,σΩ)
Ωₙ = rand(NDΩ,lenT)
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
plt_etaG = plot(t,ηᴳ,xlab="t [s]", ylab="η [m]", lw=2)

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
display(pltΩ)
savefig(pltΩ, joinpath(stats_path,"om_tilde.png"))
savefig(pltΩ, joinpath(stats_path,"om_tilde.svg"))
display(hplt_Ω)
savefig(hplt_Ω, joinpath(stats_path,"om_tilde_dist.png"))
savefig(hplt_Ω, joinpath(stats_path,"om_tilde_dist.svg"))
display(pltT2)
savefig(pltT2, joinpath(stats_path,"T2-.png"))
savefig(pltT2, joinpath(stats_path,"T2-.svg"))
display(pltΔT)
savefig(pltΔT, joinpath(stats_path,"DTev.png"))
savefig(pltΔT, joinpath(stats_path,"DTev.svg"))
display(pltFair)
savefig(pltFair, joinpath(stats_path,"FTstats.png"))
savefig(pltFair, joinpath(stats_path,"FTstats.svg"))
display(pltGEE1)
savefig(pltGEE1, joinpath(stats_path,"EWGpars3D.png"))   
savefig(pltGEE1, joinpath(stats_path,"EWGpars3D.svg"))
display(pltGEE2)
savefig(pltGEE2, joinpath(stats_path,"To-Ao.png"))
savefig(pltGEE2, joinpath(stats_path,"To-Ao.svg"))
display(hplt_Ao)
savefig(hplt_Ao, joinpath(stats_path,"Ao_dist.png"))
savefig(hplt_Ao, joinpath(stats_path,"Ao_dist.svg"))
display(hplt_To)
savefig(hplt_To, joinpath(stats_path,"To_dist.png"))
savefig(hplt_To, joinpath(stats_path,"To_dist.svg"))
display(pltGEE3)
savefig(pltGEE3, joinpath(stats_path,"tc-To.png"))
savefig(pltGEE3, joinpath(stats_path,"tc-To.svg"))
display(pltGEE4)
savefig(pltGEE4, joinpath(stats_path,"tc-Ao.png"))
savefig(pltGEE4, joinpath(stats_path,"tc-Ao.svg"))
display(pltGEE5)
savefig(pltGEE5, joinpath(stats_path,"t-beta.png"))
savefig(pltGEE5, joinpath(stats_path,"t-beta.svg"))
display(pltGEE6)
savefig(pltGEE6, joinpath(stats_path,"tc_diff.png"))
savefig(pltGEE6, joinpath(stats_path,"tc_diff.svg"))
display(hplt_Dtc)
savefig(hplt_Dtc, joinpath(stats_path,"tc_diff_dist.png"))
savefig(hplt_Dtc, joinpath(stats_path,"tc_diff_dist.svg"))
display(pltSP2)
savefig(pltSP2, joinpath(stats_path,"sp_2nd-.png"))
savefig(pltSP2, joinpath(stats_path,"sp_2nd-.svg"))
display(pltPHI2)
savefig(pltPHI2, joinpath(stats_path,"phi_2nd-.png"))
savefig(pltPHI2, joinpath(stats_path,"phi_2nd-.svg"))
display(plt_eta2)
savefig(plt_eta2, joinpath(stats_path,"eta_2nd-.png"))
savefig(plt_eta2, joinpath(stats_path,"eta_2nd-.svg"))
display(hplt_TA)
savefig(hplt_TA, joinpath(stats_path,"hist_To-Ao.png"))
savefig(hplt_TA, joinpath(stats_path,"hist_To-Ao.svg"))
display(hplt_tcA)
savefig(hplt_tcA, joinpath(stats_path,"hist_tc-Ao.png"))
savefig(hplt_tcA, joinpath(stats_path,"hist_tc-Ao.svg"))
display(hplt_tcT)
savefig(hplt_tcT, joinpath(stats_path,"hist_tc-To.png"))
savefig(hplt_tcT, joinpath(stats_path,"hist_tc-To.svg"))
display(pltG)
savefig(pltG, joinpath(stats_path,"envelopes.png"))
savefig(pltG, joinpath(stats_path,"envelopes.svg"))
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
display(plt_demo)
savefig(plt_demo, joinpath(stats_path,"carrier.png"))
savefig(plt_demo, joinpath(stats_path,"carrier.svg"))
display(plt_etaG)
savefig(plt_etaG, joinpath(stats_path,"DWE.png"))
savefig(plt_etaG, joinpath(stats_path,"DWE.svg"))
display(plt_CStats)
savefig(plt_CStats, joinpath(stats_path,"fairten_stats.png"))
savefig(plt_CStats, joinpath(stats_path,"fairten_stats.svg"))
if regress
    display(pltCONV)
    savefig(pltCONV, joinpath(stats_path,"convergence.png"))
    savefig(pltCONV, joinpath(stats_path,"convergence.svg"))
end
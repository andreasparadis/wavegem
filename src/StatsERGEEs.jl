# Load necessary modules
using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles
using FFTW, DSP, Statistics, BSplineKit
import ColorSchemes.darkrainbow

gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))
cb = darkrainbow

# Include necessary scripts for functions
include("func/text_process.jl"); include("func/signal_processing.jl")
include("func/peak_detect.jl");  include("func/runge_kutta.jl")
include("func/fugees.jl");  include("func/directories.jl")
include("func/wave_theory.jl");  include("func/jonswap.jl")

Tₚⱼ = 10.0
Hₛⱼ = 8.0
case_id = js_case(Hₛⱼ, Tₚⱼ, 3.3, true)
run_id = 4

libpath = joinpath(pwd(),"library")
case_str = "MaxFair"

pltG = plot(xlab = L"\overline{t}", ylab = L"\overline{A}", palette=:darkrainbow)
pltGEE1 = plot3d(xlab = L"\overline{T}", ylab = L"\overline{t}_c", zlab = L"\overline{A}", legend=:false)
pltGEE2 = plot(xlab = L"\overline{T}", ylab = L"\overline{A}", legend=:false)
pltGEE3 = plot(xlab = L"\overline{t}_c", ylab = L"\overline{T}", legend=:false)
pltGEE4 = plot(xlab = L"\overline{t}_c", ylab = L"\tilde{\beta}"*" [rad]", legend=:false)

t̅ₑ = 2^9/Tₚⱼ
Nₜ = 2^12
dt = t̅ₑ/(Nₜ-1)
t̅ = zeros(Float64,Nₜ)
[t̅[i] = (i-1)*dt-t̅ₑ/2 for i ∈ 1:Nₜ]

allG = zeros(Float64,Nₜ,3*11)
allΩ = zeros(Float64,3*11)
allωₚ = zeros(Float64,3*11)
count = 0

for run_id ∈ 4:6
    Decpath = joinpath(libpath,"SE","$(case_id)","$(run_id)","Decomposition")
    for evID ∈ 1:11
        evdir = joinpath(case_str,"EV$evID")
        Eventpath = joinpath(Decpath,evdir)
        GRrespath = joinpath(Eventpath,"ERGEEs")
        fid_pEWG = joinpath(GRrespath,"EWG_norm_pars")
        fid_pG = joinpath(GRrespath,"G_pars")

        # Read parameters of Gaussian Regression of wave event
        cont = parse_fxw(fid_pEWG, 1)   # Read EWGs parameters
        t̅ᶜ, T̅ₒ, A̅ₒ, β̃ = cont[:,1], cont[:,2], cont[:,3], cont[:,4]
        lenT = length(t̅ᶜ)
        cont2 = parse_fxw(fid_pG, 1)     # Read event and GR parameters
        Hₛ, Tₚ, allΩ[count*11+evID] = cont2[1:3]
        allωₚ[count*11+evID] = 2π/Tₚ

        # Resulting Gaussian approximation of the envelope
        G̅ = gauss_fun(t̅, A̅ₒ, t̅ᶜ, T̅ₒ)
        allG[:,count*11+evID] = G̅[:]

        plot!(pltG, t̅, G̅, lab=:false, line=:dot)
        plot!(pltGEE1, [T̅ₒ], [t̅ᶜ], [A̅ₒ], seriestype=:scatter)
        plot!(pltGEE2, [T̅ₒ], [A̅ₒ], seriestype=:scatter, xlim=(0,xlims(pltGEE1)[2]), ylim=(0,zlims(pltGEE1)[2]))
        plot!(pltGEE3, [t̅ᶜ], [T̅ₒ], seriestype=:scatter)
        plot!(pltGEE4, [t̅ᶜ], [β̃], seriestype=:scatter)
        # plot!(pltGEE4, [t̅ᶜ], [unwrap(β̃)], seriestype=:scatter)
    end
    global count = count+1
end

Gmean = mean!(ones(Nₜ,1),allG)
plot!(pltG,t̅, Gmean, lw=3, lab="Mean")

plt_recon = plot(t̅*Tₚⱼ, real(Hₛⱼ*Gmean.*exp.(-1im * mean(allΩ)*2π/Tₚⱼ * t̅*Tₚⱼ)))

pltΩ = plot([allΩ], seriestype=:scatter, legend=:false ,xlab="Event", ylab=L"\Omega")
# plot!(ylim=(0,2))

display(pltΩ)
display(pltGEE1)
display(pltGEE2)
display(pltGEE3)
display(pltGEE4)
display(pltG)
display(plt_recon)
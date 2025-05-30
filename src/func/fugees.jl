## Define the Elementary Gaussian model function
function elem_wg(t, a, tᶜ, T, n)
    g = a[n] * exp.(-((t.-tᶜ[n])/T[n]).^2)

    return g
end

function cgaussian(t, a, tᶜ, T)
    g = a * exp.(-((t.-tᶜ)/T).^2)

    return g
end

# Define the Gaussian model function
function gauss_fun(t, a, tᶜ, T)
    N = length(tᶜ)
    G = sum(elem_wg(t, a, tᶜ, T, n) for n in 1:N)

    return G
end

function dGdL(t, a, tᶜ, T)
    N = length(tᶜ)
    dGdL = [-2 * (t.-tᶜ[n]).^2 / T[n].^3 .* elem_wg(t, a, tᶜ, T, n) for n ∈ 1:N]

    return dGdL
end

function dLdτ(τ,T, t, A, tᶜ,u, dt)
    N = length(tᶜ)
    Gder = dGdL(t, A, tᶜ, T)
    I = [(gauss_fun(t, A, tᶜ, T) - u).*Gder[n] for n ∈ 1:N]
    dLdτ = [sum(I[n]*dt) for n ∈ 1:N]
    return dLdτ
end 

function fcsd_wg(fⱼ, Sⱼ, t, Aₒ, tᶜ, ξ)
    Nₛ = length(fⱼ)
    M = length(t)

    ωⱼ = 2π * fⱼ
    dfⱼ = abs(fⱼ[end]-fⱼ[1])/(Nₛ-1)
    η̂ = zeros(Float64,Nₛ)

    # Random.seed!(Int.(round.(rand()*1000)))
    # distr = Rayleigh(1) # Rayleigh distribution with unit scale (i.e. Rayleigh(1))
    # ξ = π/2 * rand(Uniform(-1,1))
    # α = rand(distr,Nₛ)

    # Spectrum energy through integration == Variance m₀
    m₀ = sum(Sⱼ*dfⱼ)

    ηg = zeros(Float64,M)
    for n ∈ 1:Nₛ
        η̂[n] = Aₒ/m₀ * Sⱼ[n]*dfⱼ  
        for m ∈ 1:M
        ηg[m] = ηg[m] + η̂[n] * cos(-ωⱼ[n]*(t[m].- tᶜ) + ξ)
        # ηg[m] = ηg[m] + α[n] * η̂[n] * cos(-ωⱼ[n]*(t[m].-tᶜ)+ξ)
        end
    end

    return ηg
end

function recon(prop, A̅ₒ, t̅ᶜ, T̅ₒ, β̃, Hₛ, Tₚ, Ω, t)
    ωₚ = 2π/Tₚ
    # Gaussian EWG envelope
    gₙ = cgaussian(t, A̅ₒ*Hₛ, t̅ᶜ*Tₚ, T̅ₒ*Tₚ)

    # Propagated EWG
    if prop == 0
        # Focused WG (FWG)
        ηᶜ = exp.(-1im * Ω*ωₚ * t)
        # ηᶜ,_,_ = wave_group(1,2π/ω₁,1,2π/ω₂,d,t,0, 0)
        # ηᶜ,_ = wave_group(1,Tₚ,1,2π/(Ω*ωₚ),d,t,0, 0)
        FR, MAG, ang, df,_ = one_side_asp(real(gₙ.*ηᶜ),t)
        Sᵥ = 0.5/df .* MAG.^2
        ηₙ = fcsd_wg(FR, Sᵥ, t, A̅ₒ*Hₛ, t̅ᶜ*Tₚ, β̃)
    elseif prop == 1
        # Direct Amplitude Modulation (DAM) 
        ηᶜ = exp.(-1im * Ω*ωₚ* t)
        # ηᶜ = exp.(-1im * (Ω*ωₚ .+ ωg[n]/2/π*sin.(ωg[n]*t.+π/2)) .* t)
        ηₙ = gₙ .* ηᶜ
        FR, MAG, ang, df,_ = one_side_asp(real(ηₙ) ,t)
    elseif prop == 2
        # Double Amplitude Modulation (2AM)
        ηᶜ,_,_ = wave_group(1,Tₚ,1,2π/(Ω*ωₚ),d,t,0, 0)
        ηₙ = gₙ .* ηᶜ
        FR, MAG, ang, df,_ = one_side_asp(real(ηₙ) ,t)
    else
        # Alternative 2AM
        ηᶜ,Tᵍ,Tʷ = wave_group(1,2π/ω₁,1,2π/ω₂,d,t,t̅ᶜ*Tₚ, 0)
        ηₙ = gₙ .* ηᶜ
        FR, MAG, ang, df,_ = one_side_asp(real(ηₙ) ,t)
    end

    Nₜ = length(t)

    for i ∈ 1:Nₜ
        if abs(ηₙ[i]) < 1e-6
            ηₙ[i] = 0
        end
    end    

    return gₙ, ηᶜ, ηₙ, FR, MAG
end

function regress_gauss_fun(t,u,MinPeakVal,MinPeakDist,ϵ,Nτ,dτ,Tlb,Tub)
    # Regression of a curve using a sum of Gaussian functions: G=∑gₙ, gₙ=aₙ exp(-((t-tᶜₙ)/Tₙ)²)
    # Required inputs:
    #   u=u(t), t: vectors of independent and dependent variables 
    #   MinPeakVal, MinPeakDist: thresholds for peak detection
    #   ϵ,Nτ,dτ,Tlb,Tub: Runge-Kutta solution parameters = accuracy, max. No of iterations, step of fictitious time, constraints for T

    # Initialise parameter vectors
    Ø = Array{Float64}(undef,0)
    Tₒ, Aₒ, tᶜₒ = [Ø[:] for _ = 1:3]
    i⁺ₒ = Array{Int64}(undef,0)

    # Peak detection of maxima above given threshold and distance between them
    # the "true" argument sorts in descending order based on A
    A, tᶜ, i⁺, uₜ, uₜₜ = peaks_max_ext(u,t, MinPeakVal, MinPeakDist,true)
    N = length(tᶜ)

    ## Initial condition for L based on ∂ₜ²g(tᶜ) = ∂ₜ²u(tᶜ)
    T = sqrt.(abs.(2*A./uₜₜ[i⁺]))     

    # Definition of the ODE's dLₘ 
    τ = zeros(Float64,Nτ)
    [τ[i] = (i-1)*dτ for i ∈ 2:Nτ]
    println("τₑ = ", τ[end])
    fun(τ,T) = dLdτ(τ,T, t, A, tᶜ,u, dt)

    # Solve optimization problem
    Tₒₚₜ = RK4_nln_sys(τ,T,fun,Tlb,Tub,ϵ)

    # Optimal values
    pltCONV = plot(palette=:darkrainbow)

    for n ∈ 1:N
        plot!(pltCONV, Tₒₚₜ[2][:,n],xlab="i",ylab=L"T_o^i [s]",xscale=:log10,lab="To[$(n)]",lw=2,leg=:outerbottom,legendcolumns=N)
        if Tₒₚₜ[1][n] > Tcut
            push!(Tₒ, Tₒₚₜ[1][n])
            push!(Aₒ, A[n])
            push!(tᶜₒ, tᶜ[n])
            push!(i⁺ₒ,i⁺[n])
        end
    end

    # Resulting Gaussian approximation of the envelope
    G = gauss_fun(t, Aₒ, tᶜₒ, Tₒ)

    return Tₒ, Aₒ, tᶜₒ, i⁺ₒ, G, pltCONV
end
function dispersion(T,d)
    g::Float64 = 9.81
    Ø = Array{Float64}(undef,0)
    λ = Ø[:]

    λ₀ = g/(2*π) * T^2
    push!(λ, λ₀) 

    i = 1
    ϵ = 1.0
    while ϵ>1e-3 && i<100
        push!(λ, λ₀ * tanh(2*π/λ[i]*d))
        ϵ = abs((λ[i+1]-λ[i])/λ[i])
        i = i + 1
    end

    if i == 1000
        print("The value of λ has not converged after a 1000 iterations. The final relative error is ", ϵ, "\n")
    end

    λₑ = λ[end]

    return λₑ, λ, ϵ, i
end

function wave_qnts(T,d)
    f = 1/T
    ω = 2π/T
    λ₀ = g/(2*π) * T^2 # Deep water wavelength / Initial guess
    
    if d > λ₀/2
        λ = λ₀
        println("Deep water (T=$(round(T*100)/100) [s], λ=$(round(λ*100)/100) [m])")
    elseif d>λ₀/20 && d<=λ₀/2 
        λ,_ = dispersion(T,d)
        println("Intermediate water (T=$(round(T*100)/100) [s], λ=$(round(λ*100)/100) [m])")
    else
        λ,_ = dispersion(T,d)
        println("Shallow water (T=$(round(T*100)/100) [s], λ=$(round(λ*100)/100) [m])")
    end

    κ = 2π/λ
    υᶜ = ω/κ

    wave = Wave(f, T, ω, λ, κ, υᶜ)

    return wave
end

function rand_sea(fⱼ, Sⱼ, Tₑ, dt, Nₜ, phi_id, A_id, Ldom, d, Cr, runstr)
    Nₛ = length(fⱼ)

    ω = 2π * fⱼ;    k = 1/g * ω.^2;     #υ = ω./k; υ[1]=0

    sp_peak = findmax(Sⱼ);     iₚ = sp_peak[2];    fₚ = fⱼ[iₚ]

    for i ∈ 1:iₚ
        if 2π/ω[i] < Tₑ
            λ,_ = dispersion(2π/ω[i],d)
            k[i] = 2π/λ
        end
    end
    υ = sqrt.(g./k .* tanh.(k*d)); 

    dfⱼ = abs(fⱼ[end]-fⱼ[1])/(Nₛ-1)

    # Spatial discretisation
    ## x-axis
    dx = round(dt*υₚ/Cr *10)/10
    Cr = dt*υₚ/dx
    Nₓ = Int64(round(Ldom/dx))+1
    ## z-axis   
    dz = 2*dx
    Nz = Int64(round(d/dz))+1

    # Initialize vectors
    t = zeros(Float32,Nₜ)
    [t[i] = (i-1)*dt for i ∈ 1:Nₜ]

    x = zeros(Float64,Nₓ)
    [x[i] = (i-1)*dx for i ∈ 1:Nₓ]

    z = zeros(Float64,Nz)
    [z[i] = -(i-1)*dz for i ∈ 1:Nz]

    𝚶⃗ₜ = zeros(Float32,Nₜ)
    𝚶⃗ₛ = zeros(Float32,Nₛ)
    𝚶⃗ₓ = zeros(Float32,Nₓ)

    ηₜ, ηₜᵈ, u, u̇, ẇ, Φ, U = (𝚶⃗ₜ[:] for _ = 1:9)
    U̅ = 𝚶⃗ₛ[:]
    ηₓ = 𝚶⃗ₓ[:]

    # Wave component amplitudes
    η̂ = sqrt.(2*Sⱼ*dfⱼ) 
    m₀ = Hₛ^2/16
    H₀ = sqrt(2*m₀*log(Nₛ))

    # Phase
    amp = rand(Normal(),Nₛ) + 1im*rand(Normal(),Nₛ) # Random complex vector (Normal distribution, 0 mean, 1 variance)
    if phi_id == 0
        ϕ = angle.(amp) # ϕ ∈ [0,2π]
    else
        cont = parse_fxw_pf(runstr*"/$(phi_id-1)/phi", 0, 1)
        ϕ = cont[:,2] .- π/2
        print("-$(phi_id)π/2 phase shift \n")
    end

    # Randomly distributed amplitude coefficient
    if A_id == 0
        A = ones(Float64,Nₛ)
        println("Unit amplitude coefficients")
    elseif A_id == 1
        A = rand(Normal(),Nₛ)
        println("Normally distributed amplitude coefficients")
    elseif A_id == 2
        seed = Int64(round(rand(1)[1]*1000));   Random.seed!(seed);    
        distr = Rayleigh(1);    A = rand(distr,Nₛ) # Rayleigh distribution with unit scale
        println("Rayleigh distributed amplitude coefficients")
    else
        A = abs.(real(amp))
        println("Amplitude coefficients from normally distributed real part of amp")
    end

    # Surface Wave kinematics at x=0
    Ĥ = A .* η̂
    Û = ω .* Ĥ
    α̂ = ω.^2 .* Ĥ
    Φ̂ = g * Ĥ ./ ω

    for m ∈ 1:Nₛ
        θ = -ω[m]*t .+ ϕ[m]
        ϵ = k[m]*x .+ ϕ[m]
        csθ = cos.(θ)
        snθ = sin.(θ)
        csϵ = cos.(ϵ)

        ηₜ[:] += Ĥ[m] * csθ
        ηₜᵈ[:] += Ĥ[m] * exp(-k[m]*d) * csθ
        Φ[:] += Φ̂[m] *snθ
        u[:] += Û[m] * csθ
        u̇[:] += α̂[m] * snθ
        ẇ[:] += - α̂[m] * snθ

        uᶻ = ω[m] * η̂[m] * exp.(k[m]*z)
        U̅[m] = mean(uᶻ)
        U[:] += U̅[m] * csθ

        ηₓ[:] += Ĥ[m] * csϵ
    end

    return η̂, t, dt, ηₜ, u, u̇, ẇ, Φ, ηₜᵈ, ϕ, x, dx, Nₓ, ηₓ, z, dz, Nz, U
end
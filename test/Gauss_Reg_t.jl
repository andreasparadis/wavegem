## Define the Elementary Gaussian model function
function elem_wg(t, a, tᶜ, T, n)
    g = a[n] * exp.(-((t.-tᶜ[n])/T[n]).^2)

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
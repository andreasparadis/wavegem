function spectrum(Hₛ,Tₚ,γ,Tᵢ,Tₑ,N)
    g::Float64 = 9.81
    S = zeros(Float64,N)

    if Tᵢ<=0
        print("\n -------------------- \n
WARNING: The lowest wave period must be larger than 0.\n
The value of Tᵢ was changed to Tᵢ=0.4 s\n")
        Tᵢ = Float64(Tᵢ)
        Tᵢ = 0.4
    end

    fₑ = 1/Tₑ;    ωᵢ = 2π/Tₑ
    fᵢ = 1/Tᵢ;    ωₑ = 2π/Tᵢ
    f = range(fₑ,fᵢ,N)
    ω = 2π*f

    ωₚ = 2π/Tₚ # Peak circular frequency
    fₚ = 1/Tₚ # Peak frequency
    υₚ = g/ωₚ # wave celerity at deep water

    α = ((Hₛ*ωₚ^2)/(2.227g))^2
    x̃ = (ωₚ/22)^-3
    U₁₉ = υₚ/1.14
    U₁₀ = U₁₉/1.026
    x = x̃*U₁₀^2/g

    m₀ = Hₛ^2/16 # m₀ = 0.31*α*g^2 / ωₚ^4
    m₋₁ = 5.655*ωₚ^-1 * m₀
    m₁ = 0.191*ωₚ * m₀
    m₂ = 0.043*ωₚ^2 * m₀
    m₃ = 0.013*ωₚ^3 * m₀

    T1 = 2π*m₀/m₁
    Tz = 2π*sqrt(m₀/m₂)

    A = α*g^2 / (2π)^4
    B = 5/4*fₚ^4

    for i ∈ 1:N
        if ω[i] < ωₚ
            σ = 0.07
            r = exp(-1/2*((f[i]-fₚ)/(σ * fₚ))^2)
        else
            σ = 0.09
            r = exp(-1/2*((f[i]-fₚ)/(σ * fₚ))^2)
        end

        S[i] = A/f[i]^5 * exp(-B/f[i]^4) * γ^r
    end

    # JS_pars = []

    println("----------------------------")
    println("JONSWAP spectrum parameters:")
    println("Tₚ=",Tₚ," [s]")
    println("Hₛ=",Hₛ," [m]")
    println("γ=",γ," [-]")
    println("α=",α," [-]")
    # println("x=",x," [m]")
    # println("U₁₀=",U₁₀," [m/s]")
    println("m₀=",m₀," [m²]")
    # println("υₚ=",υₚ," [m/s]")
    println("T₁=",T1,"[s]")
    println("Tz=",Tz,"[s]")
    println("----------------------------")

    return f,S
end

function height_scatter(Tᵢ,Tₑ,N)
    Tₚ = range(Tᵢ,Tₑ,N)
    Hₛ = 0.101*Tₚ.^1.67

    return Tₚ,Hₛ
end
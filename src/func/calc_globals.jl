function time_lims(flong, Tₑ)
    fcut = 1/Tₑ+1/2
    Tᵢ = 1/fcut             # Cut-off period of JONSWAP spectrum [s]
    λ⁻ = g/(2*π) * Tᵢ^2     # Shortest wavelength (Deep water) [m]
    υ⁻ = λ⁻/Tᵢ              # Slowest wave [m/s]

    if flong 
        tₑ = 2^2 * nextpow(2,2*Ldom/υ⁻) # Duration of long simulation [s]
    else
        tₑ = nextpow(2,2*Ldom/υ⁻) # Duration of simulation [s]
    end

    return Tᵢ, λ⁻, υ⁻, tₑ, fcut
end
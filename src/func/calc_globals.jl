function time_lims(fcut, flong)
    Tᵢ = 1/fcut             # Cut-off period of JONSWAP spectrum [s]
    λ⁻ = g/(2*π) * Tᵢ^2     # Shortest wavelength (Deep water) [m]
    υ⁻ = λ⁻/Tᵢ              # Slowest wave [m/s]

    if flong 
        tₑ = 180*60 # 3-hour duration of simulation [s]
    else
        tₑ = round(2*Ldom/υ⁻ /10)*10 # Duration of simulation [s]
    end

    return Tᵢ, λ⁻, υ⁻, tₑ
end
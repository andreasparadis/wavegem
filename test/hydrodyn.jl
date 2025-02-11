function float_stabil(T,r₁,r₂,h₂)
# Floater consisting on cylindrical body (r₁) and heave plate (r₂,h₂)
# T: Draught

    # Centre of bouyancy
    h₁ = T - h₂
    A₁ = 2*r₁*h₁
    A₂ = 2*r₂*h₂

    zᶜ₁ = -h₁/2
    zᶜ₂ = -(h₁+h₂/2)

    zᵦ = (zᶜ₁*A₁ + zᶜ₂*A₂)/(A₁+A₂)

    # Mass of body
    Vᶠ = π*r₁^2*h₁ + π*r₂^2*h₂  # Displaced fluid volume
    ρᶠ = 1000
    mᵇ = ρᶠ*Vᶠ

    # Metacenter
    Iy = π/4*r₁^4

    zₘ = zᵦ + Iy/Vᶠ

    return zᵦ, zₘ
end
function cartesian_to_spherical(xₚ)
    # xₚ: Time dependent position vector of a specific point (trajectory)
    # xₚ = [3xN] {x̂;ŷ;ẑ}, N: number of time steps
    x, y, z = xₚ[1,:], xₚ[2,:], xₚ[3,:]
    
    # Compute the radial distance r
    r = sqrt.(x.^2 + y.^2 + z.^2)
    # Compute the polar angle θ (in radians), θ is in [0, π]
    θ = acos.(z ./ r)
    # Compute the azimuthal angle φ (in radians), φ is in [-π, π], handling all quadrants
    φ = atan.(y, x)
    
    rθφ = [r';θ';φ']

    # Check
    r̂ = [sin.(θ).*cos.(φ) sin.(θ).*sin.(φ) cos.(θ)]
    rₚ = r .* r̂

    return rθφ, rₚ
end

function spherical_to_cartesian(rθφ)
    r, θ, φ = rθφ[1], rθφ[2], rθφ[3]

    r̂ = [sin(θ)*cos(φ); sin(θ)*sin(φ); cos(θ)]
    θ̂ = [cos(θ)*cos(φ); cos(θ)*sin(φ); -sin(θ)]
    ϕ̂ = [-sin(φ); cos(φ); 0]

    T = [r̂'; θ̂'; ϕ̂']

    xₚ = T * rₚ

    return xₚ
end

function derivatives(x,t)
    # Initialize variables
    m = length(t)
    𝚶⃗ₘ = zeros(Float64,m)
    ẋ, ẍ = 𝚶⃗ₘ[:], 𝚶⃗ₘ[:]

    # 1st differentiation
    for i ∈ 2:m-1
        ẋ[i] = (x[i+1] - x[i-1]) / (t[i+1]-t[i-1])
    end

    # Linear interpolation for first element
    ẋ[1] = 2*ẋ[2]-ẋ[3]

    # 2nd differentiation and storing of maxima and minima
    for i ∈ 2:m-1
        ẍ[i] = (ẋ[i+1] - ẋ[i-1]) / (t[i+1]-t[i-1])
    end

    # Linear interpolation of last element
    ẍ[end] = 2*ẍ[end-1]-ẍ[end-2]

    return ẋ, ẍ
end

function time_integration(t::Vector, ẋ::Vector, x0::Float64)
    n = length(t)
    dt = (t[end] - t[1])/(n-1)
    x = zeros(Float64, n)  # Array to store values
    x[1] = x0             # Set initial value

    # Simpson's rule for most of the points
    for i in 2:n-2
        x[i+1] = x[i-1] + (dt / 3) * (ẋ[i-1] + 4*ẋ[i] + ẋ[i+1])
    end

    return x
end

function trapez_integral(t::Vector, ẋ::Vector, x0::Float64)
    n = length(t)
    dt = (t[end] - t[1])/(n-1)
    x = zeros(Float64, n)  # Array to store values
    x[1] = x0             # Set initial value

    # Simpson's rule for most of the points
    for i in 1:n-1
        x[i+1] = x[i] + dt  * (ẋ[i] + ẋ[i+1])/2
    end

    return x
end
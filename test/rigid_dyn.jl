function kinematics_FOWT(x,y,z,t)
    # The function returns the kinematics expressed by the x̂, ŷ, ẑ components
    # x,y,z = [Nx1]: Time dependent vectors, N number of time steps 

    # Definition of angles (positive rotation using the right hand rule)
    # ϕ : roll (rotation around x-axis)
    # θ : pitch (rotation around y-axis)
    # ψ : yaw (rotation around z-axis)
    r = sqrt.(x.^2 + y.^2 + z.^2)
    ϕ = atan.(z,y) 
    θ = atan.(x,z) 
    ψ = atan.(y,x)

    # Velocities & accelerations
    u, u̇ = derivatives(x,t)
    v, v̇ = derivatives(y,t)
    w, ẇ = derivatives(z,t)

    # Angular velocities & accelerations
    ωˣ, ω̇ˣ = derivatives(unwrap(ϕ),t)
    ωʸ, ω̇ʸ = derivatives(unwrap(θ),t)
    ωᶻ, ω̇ᶻ = derivatives(unwrap(ψ),t)

    υ = [u v w]
    a = [u̇ v̇ ẇ]
    ω = [ωˣ ωʸ ωᶻ]
    ω̇ = [ω̇ˣ ω̇ʸ ω̇ᶻ]

    return r, ϕ, θ, ψ, υ, a, ω, ω̇
end

function kinematics(rₚ,t)
    # The function returns the kinematics expressed by the x̂, ŷ, ẑ components
    # rₚ: Time dependent position vector of a specific point (trajectory) 
    # rₚ = [3xN] {x̂;ŷ;ẑ}, N: number of time steps
    x, y, z = rₚ[1,:], rₚ[2,:], rₚ[3,:]
    N = length(t)

    # Spherical coordinates - Physics convention 
    # *The definition of θ and ϕ are not necessarily the same as in the caller script,
    # *take care when using the returned values or only use within the function
    r = sqrt.(x.^2 + y.^2 + z.^2)
    θ = acos.(z ./ r) 
    ϕ = atan.(y, x) 
    
    # Unit vectors
    r̂ = [sin.(θ).*cos.(ϕ) sin.(θ).*sin.(ϕ) cos.(θ)]
    θ̂ = [cos.(θ).*cos.(ϕ) cos.(θ).*sin.(ϕ) -sin.(θ)]
    ϕ̂ = [-sin.(ϕ) cos.(ϕ) zeros(Float64,N)]

    # Calculate time derivatives
    ṙ, r̈ = derivatives(r,t)
    θ̇, θ̈ = derivatives(θ,t)
    ϕ̇, ϕ̈ = derivatives(ϕ,t)

    υ = zeros(Float64,3,N)
    a = zeros(Float64,3,N)
    ω = zeros(Float64,3,N)
    Φ = zeros(Float64,3,N)
    Θ = zeros(Float64,3,N)
    Ψ = zeros(Float64,3,N)

    for i ∈ 1:N
        T = [r̂[N,:] θ̂[N,:] ϕ̂[N,:]]
        x̂ = T[1,:]' 
        ŷ = T[2,:]' 
        ẑ = T[3,:]' 
        # Velocity
        υ[:,i] = [ṙ[i] r[i]*θ̇[i] r[i]*sin(θ[i]).*ϕ̇[i]] * T'

        # Acceleration
        aᵣ = r̈[i].-r[i]*θ̇[i]^2 - r[i]*ϕ̇[i]^2 *sin(θ[i])^2
        aₜ = r[i]*θ̈[i]+2*ṙ[i]*θ̇[i]-r[i]*ϕ̇[i]^2 *sin(θ[i])*cos(θ[i])
        aₚ = r[i]*ϕ̈[i]*sin(θ[i])+2*ṙ[i]*ϕ̇[i]*sin(θ[i])+2*r[i]*θ̇[i]*ϕ̇[i]*cos(θ[i])
        a[:,i] = [aᵣ aₜ aₚ] * T'

        # Angular velocities
        ωˣ = θ̇[i]*cos(ϕ[i]) + ϕ̇[i]*sin(θ[i])*sin(ϕ[i])
        ωʸ = θ̇[i]*sin(ϕ[i]) - ϕ̇[i]*sin(θ[i])*cos(ϕ[i])
        ωᶻ = ϕ̇[i]
        ω[:,i] = [ωˣ ωʸ ωᶻ]
    end

    return υ, a, ω, r, θ, ϕ
end

function findCG(zG)
    rem = zeros(Float64,3,Nₜ)
    zG0 = range(163,164,100)
    smin = zeros(Float64,100)

    for j ∈ 1:100
        xG0 = zG0[j]*sin(θ₀)
        yG0 = zG0[j]*sin(ϕ₀)

        # Center of gravity motion
        xG = xG0 .+ dx
        yG = yG0 .+ dy
        zG = zG0[j] .+ dz
        
        for i ∈ 1:Nₜ
            Ω = [ωˣₐ[i]; ωʸₐ[i]; ωᶻₐ[i]]
            rᴳ = [X[i,2]; Y[i,2]; Z[i,2]] .- [xG[i]; yG[i]; zG[i]]

            rem[:,i] = cross(Ω,rᴳ) .+ [u₀[i]; v₀[i]; w₀[i]] .- [I.u[i,2]; I.v[i,2]; I.w[i,2]]
        end
        smin[j] = sum(abs.(rem[1,:]))
    end

    plot(zG0,smin)
    plot(t, rem[1,:])

    Val, Pos, _ = peaks_extend(smin,zG0, 0, 0)

    # # Alternatively
    # Mark = 2
    # uᴵ₁ = I.u[:,Mark]
    # x₁ = -Z[:,2].*sin.(B.θ[:,Mark])
    # y₁ = zeros(Float64,Nₜ)
    # z₁ = (uᴵ₁ .- u₀ .+ B.ωᶻ[:,Mark].*y₁)./B.ωʸ[:,Mark]
    # R₁ = sqrt.(x₁.^2 .+ y₁.^2 .+ z₁.^2)

    # FOWT
    # r₀ = [0;0;0]
    # if fdate == "2024-10-15"
    #     r₀ = [168.451; -937.256; -169.937] # 2024-10-15
    # elseif fdate == "2024-10-14"
    #     r₀ = [167.357; -938.395; -168.659] # 2024-10-14
    # elseif fdate == "2024-10-10"
    #     r₀ = [165.879; -932.311; -163.647] # 2024-10-10
    # elseif fdate == "2024-10-09"
    #     r₀ = [165.019; -933.612; -160.699] # 2024-10-09
    # elseif fdate == "2024-10-08"
    #     r₀ = [165.673; -933.149; -159.757] # 2024-10-08
    # else
    #     r₀ = [168.752; -937.946; -136.84] # 2024-10-07
    # end
end


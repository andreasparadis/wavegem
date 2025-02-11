mutable struct pos 
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

mutable struct Markers 
    M1::pos
    M2::pos
    M3::pos
    M4::pos
    M5::pos
    M6::pos
    M7::pos
    M8::pos
    M9::pos
    M10::pos
end

mutable struct kin 
    x::Matrix{Float64}
    y::Matrix{Float64}
    z::Matrix{Float64}
    ϕ::Matrix{Float64}
    θ::Matrix{Float64}
    ψ::Matrix{Float64}
    u::Matrix{Float64}
    v::Matrix{Float64}
    w::Matrix{Float64}
    ωˣ::Matrix{Float64}
    ωʸ::Matrix{Float64}
    ωᶻ::Matrix{Float64}
    u̇::Matrix{Float64}
    v̇::Matrix{Float64}
    ẇ::Matrix{Float64}
    ω̇ˣ::Matrix{Float64}
    ω̇ʸ::Matrix{Float64}
    ω̇ᶻ::Matrix{Float64}
end

mutable struct point_kin 
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    ϕ::Vector{Float64}
    θ::Vector{Float64}
    ψ::Vector{Float64}
    u::Vector{Float64}
    v::Vector{Float64}
    w::Vector{Float64}
    ωˣ::Vector{Float64}
    ωʸ::Vector{Float64}
    ωᶻ::Vector{Float64}
    u̇::Vector{Float64}
    v̇::Vector{Float64}
    ẇ::Vector{Float64}
    ω̇ˣ::Vector{Float64}
    ω̇ʸ::Vector{Float64}
    ω̇ᶻ::Vector{Float64}
end

mutable struct Acc 
    u̇::Vector{Float64}
    v̇::Vector{Float64}
    ẇ::Vector{Float64}
    ωˣ::Vector{Float64}
    ωʸ::Vector{Float64}
    ωᶻ::Vector{Float64}
    ϕ::Vector{Float64}
    θ::Vector{Float64}
    ψ::Vector{Float64}
end
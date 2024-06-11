using Plots
using DataFrames
using LaTeXStrings
using LinearAlgebra

g::Float64 = 9.81
n::Int64 = 1000
T::Float64 = 2
H::Float64 = 2
η = zeros(Float64,n)
η₁ = zeros(Float64,n)
η₂ = zeros(Float64,n)
η₃ = zeros(Float64,n)

ω₁ = 2*π/T
κ₁ = ω₁^2/g
λ₁ = 2*π/κ₁

ω₂ = 1.1*ω₁
κ₂ = ω₂^2/g
λ₂ = 2*π/κ₂

ω₃ = 1.2*ω₁
κ₃ = ω₃^2/g
λ₃ = 2*π/κ₃

η̂ = H/2
κ₊ = (κ₁+κ₂)/2
ω₊ = (ω₁+ω₂)/2
λ₊ = 2*π/κ₊
T₊ = 2*π/ω₊

κ₋ = abs(κ₁-κ₂)/2
ω₋ = abs(ω₁-ω₂)/2
λ₋ = 2*π/κ₋
T₋ = 2*π/ω₋


x = range(-10*λ₋,10*λ₋,n)
t = range(0,5*T₋,n)

for i=1:n
    # η₁[i] = η̂ * cos(κ₁*x[i])
    η₁[i] = η̂ * cos(-ω₁*t[i])

    # η₂[i] = η̂ * cos(κ₂*x[i])
    η₂[i] = η̂ * cos(-ω₂*t[i])

    # η₃[i] = η̂ * cos(κ₃*x[i])
    η₃[i] = η̂ * cos(-ω₃*t[i])

    η[i] = η₁[i] + η₂[i] #+ η₃[i]
end

# plot(x./λ₊, η₁, xlabel = L"x/λ", ylabel = L"\eta", label = L"\eta_1", linestyle=:dot)
# plot!(x./λ₊, η₂, xlabel = L"x/λ", ylabel = L"\eta", label = L"\eta_2", linestyle=:dot)
# plot!(x./λ₊, η, xlabel = L"x/λ", ylabel = L"\eta", label = L"\eta", linewidth=3)

plot(t./T₋, η₁, xlabel = "t/T", ylabel = L"\eta", label = L"\eta_1")
plot!(t./T₋, η₂, xlabel = "t/T", ylabel = L"\eta", label = L"\eta_2")
plot!(t./T₋, η, xlabel = "t/T", ylabel = L"\eta", label = L"\eta", linewidth=3)
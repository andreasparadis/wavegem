using Plots
using DataFrames

n::Int16 = 1000
a::Float32 = 10.0
x0::Float32 = 1.0

ℯ=exp(1)

x = range(0,4,n)
f1 = zeros(Float32,n,1)
f2 = zeros(Float32,n,1)

for i = 1:n
    f1[i] = ℯ^(-a/x[i]^4)
    f2[i] = 0.5*(16*a^2/x0^10 - 20*a/x0^6) * ℯ^(-a/x0^4) * (x[i]-x0)^2 + 4*a/x0*ℯ^(-a/x0^4)*(x[i]-x0) + ℯ^(-a/x0^4)
end

plt = plot(x, f1, xlabel = "x", ylabel = "f(x)", 
       label = "Original")

plt2 = plot!(x, f2, xlabel = "x", ylabel = "f(x)", 
       label = "Taylor expansion")

display(plt)
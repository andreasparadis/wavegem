# Functions implementing various implementations of Runge-Kutta methods of solution for ODEs
function RK4(x,y0,f, scheme)

    L = length(x)
    Δx = (x[end]-x[1])/(L-1)
    y = zeros(Float64,L)
    y[1] = y0

    ## Choose a specific scheme for the RK4 implementation
    if scheme == "Gill"
        ### Weight coefficients of specific scheme
        w = [1/6; 1/3*(1-1/sqrt(2)); 1/3*(1+1/sqrt(2)); 1/6]

        for i ∈ 1:L-1
            k₁ = Δx * f(x[i], y[i])
            k₂ = Δx * f(x[i]+Δx/2, y[i]+k₁)
            k₃ = Δx * f(x[i]+Δx/2, y[i] + (-1/2+1/sqrt(2))*k₁ + (1-1/sqrt(2))*k₂)
            k₄ = Δx * f(x[i+1], y[i] - 1/sqrt(2)*k₂ + (1+1/sqrt(2))*k₃)

            ### Calculate the value of the next step
            y[i+1] = y[i] + w[1]*k₁ + w[2]*k₂ + w[3]*k₃ + w[4]*k₄
        end
    else
        w = [1/6; 1/3; 1/3; 1/6]

        for i ∈ 1:L-1
            k₁ = Δx * f(x[i], y[i])
            k₂ = Δx * f(x[i]+Δx/2, y[i]+k₁/2)
            k₃ = Δx * f(x[i]+Δx/2, y[i]+k₂/2)
            k₄ = Δx * f(x[i+1], y[i]+k₃)

            y[i+1] = y[i] + w[1]*k₁ + w[2]*k₂ + w[3]*k₃ + w[4]*k₄
        end
    end

    return y
end

function RK4_nln_sys(x,y0,f)

    L = length(x)
    W = length(y0)
    Δx = (x[end]-x[1])/(L-1)
    y = zeros(Float64,L,W)
    y[1,:] = y0

    w = [1/6; 1/3; 1/3; 1/6]

    for i ∈ 1:L-1
        k₁ = Δx * f(x[i], y[i,:])
        k₂ = Δx * f(x[i]+Δx/2, y[i,:].+k₁/2)
        k₃ = Δx * f(x[i]+Δx/2, y[i,:].+k₂/2)
        k₄ = Δx * f(x[i+1], y[i,:].+k₃)

        y[i+1,:] = y[i,:] + w[1]*k₁ + w[2]*k₂ + w[3]*k₃ + w[4]*k₄
    end

    return y[end,:], y
end
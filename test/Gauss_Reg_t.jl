## Define the Elementary Gaussian model function
function elem_wg(t, a, tᶜ, T, n)
    g = a[n] * exp.(-((t.-tᶜ[n])/T[n]).^2)

    return g
end

# Define the Gaussian model function
function gauss_fun(t, a, tᶜ, T)
    N = length(tᶜ)
    G = sum(elem_wg(t, a, tᶜ, T, n) for n in 1:N)

    return G
end

function dGdL(t, a, tᶜ, T)
    N = length(tᶜ)
    dGdL = [-2 * (t.-tᶜ[n]).^2 / T[n].^3 .* elem_wg(t, a, tᶜ, T, n) for n ∈ 1:N]

    return dGdL
end

function dLdτ(τ,T, t, A, tᶜ,u, dt)

    N = length(tᶜ)
    Gder = dGdL(t, A, tᶜ, T)
    I = [(gauss_fun(t, A, tᶜ, T) - u).*Gder[n] for n ∈ 1:N]
    dLdτ = [sum(I[n]*dt) for n ∈ 1:N]
    return dLdτ
end 
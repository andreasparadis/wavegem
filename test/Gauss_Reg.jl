## Define the Elementary Gaussian model function
function elem_wg(x, a, xᶜ, L, n)
    g = a[n] * exp.(-((x.-xᶜ[n])/L[n]).^2)

    return g
end

# Define the Gaussian model function
function gauss_fun(x, a, xᶜ, L)
    N = length(xᶜ)
    G = sum(elem_wg(x, a, xᶜ, L, n) for n in 1:N)

    return G
end

function dGdL(x, a, xᶜ, L)
    N = length(xᶜ)
    dGdL = [-2 * (x.-xᶜ[n]).^2 / L[n].^3 .* elem_wg(x, a, xᶜ, L, n) for n ∈ 1:N]

    return dGdL
end

function dLdτ(τ,L)
    global x, A, xᶜ,u, dx

    N = length(xᶜ)
    Gder = dGdL(x, A, xᶜ, L)
    I = [(gauss_fun(x, A, xᶜ, L) - u).*Gder[n] for n ∈ 1:N]
    dLdτ = [sum(I[n]*dx) for n ∈ 1:N]
    return dLdτ
end 
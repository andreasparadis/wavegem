using LsqFit
using Plots

# Define the Elementary Gaussian model function
function elementary_gaussian_model(x, p, i)
    g = 1/(p[(i-1)*3+2])*sqrt(2π) .* exp.(-((x .- p[(i-1)*3+1]) / (sqrt(2)*p[(i-1)*3+2])).^2)

    return g
end

# Define the Gaussian model function
function gaussian_model(x, p)
    n_peaks = length(p) ÷ 3
    G = sum(elementary_gaussian_model(x, p, i) for i in 1:n_peaks)

    return G
end

# Define a function to generate synthetic data with multiple peaks
function generate_multi_peak_data(x, true_parameters)
    n_peaks = length(true_parameters) ÷ 3
    return sum(gaussian_model(x, true_parameters[(i-1)*3+1:i*3]) for i in 1:n_peaks) + 0.1 * randn(size(x))
end

# Define the objective function for least squares fitting
function multi_peak_objective_function(p, x, y_observed)
    J = gaussian_model(x, p) - y_observed
    return J
end

# Generate synthetic data with an arbitrary number of peaks
n_peaks_generated = 9  # Change this to the desired number of peaks
true_parameters_multi_peak = rand(27)  # μ, σ, A for each peak

# Generate synthetic data with multiple peaks
x_data = collect(range(0, stop=12, length=1000))
y_true_multi_peak = generate_multi_peak_data(x_data, true_parameters_multi_peak)

# Perform the least squares fitting
initial_guess_multi_peak = rand(3 * n_peaks_generated)  # Initial guess for parameters
fit_result_multi_peak = curve_fit((p, x) -> multi_peak_objective_function(true_parameters_multi_peak, x_data, y_true_multi_peak) , x_data, y_true_multi_peak, initial_guess_multi_peak)
best_fit_parameters_multi_peak = fit_result_multi_peak.param
yfit = generate_multi_peak_data(x_data, best_fit_parameters_multi_peak)

# Plot the results
plot(x_data, y_true_multi_peak, label="True Data")
plot!(x_data, yfit, label="Best Fit", linewidth=2)

# # Plot individual peaks for better visualization
# for i in 1:n_peaks_generated
#     display(plot!(x_data, gaussian_model(x_data, best_fit_parameters_multi_peak[(i-1)*3+1:i*3]), label="Peak $i", linestyle=:dash))
# end
function dist_fit(dataset, dist_type, bins)
    # Sample Statistics
    x̅ = mean(dataset)   # Sample mean
    s = std(dataset)    # Sample standard deviation
    samp_params = [x̅;s]

    # Fit the specified distribution to the dataset
    fitted_dist = fit(dist_type, dataset)

    # Extract the parameters of the fitted distribution
    params = [getfield(fitted_dist, field) for field in fieldnames(typeof(fitted_dist))]

    # Create a range of values for plotting the PDF
    if dist_type == Normal
        value_range = range(x̅-4*s, x̅+4*s, length=bins)
    elseif dist_type == Weibull || dist_type == LogNormal
        value_range = range(0, maximum(dataset), length=bins)
    else
        value_range = range(minimum(dataset), maximum(dataset), length=bins)
    end

    # Calculate the PDF of the fitted distribution
    pdf_values = pdf.(fitted_dist, value_range)

    # Create a histogram of the dataset
    hst_plt = histogram(dataset, normalize=:pdf, lab=:false, palette=[cb[8]; cb[11]; cb[5]])

    # Plot the PDF of the fitted distribution
    if dist_type == Weibull
        param_lab = "$dist_type : κ=$(round(params[1]*1e3)/1e3), λ=$(round(params[2]*1e3)/1e3)"
    else
        param_lab = "$dist_type : μ=$(round(params[1]*1e3)/1e3), σ=$(round(params[2]*1e3)/1e3)"
    end

    plot!(hst_plt, value_range, pdf_values, lw=2, label=param_lab, xlab="Value", ylab="Density")

    ## Q-Q plot
    data_quantiles = sort(dataset)
    theor_quantiles = quantile.(fitted_dist, (1:length(dataset)) / (length(dataset) + 1))
    log_like = sum(log.(pdf.(fitted_dist, dataset)))

    QQplt = plot(xlabel="Theoretical Quantiles", ylabel="Sample Quantiles", title="Q-Q Plot", palette=[cb[11]; cb[8]])
    plot!(QQplt, theor_quantiles, data_quantiles, lab="Log-likelihood: $(round(log_like))", lw=2)
    plot!(QQplt, theor_quantiles, theor_quantiles, linestyle=:dash, lab="45°")

    return pdf_values, fitted_dist, params, samp_params, hst_plt, QQplt, value_range
end

function copula2d_fit_eval(CopType,X₁,X₂,f₁,f₂,res,grid_size)
    # X₁, X₂:     original samples
    # f₁, f₂:     pdfs of marginal distributions
    # res:        resolution of sampling from pdfs
    # grid_size:  resolution of the grid for the evaluation of the copula
    
    # Calculate the CDFs of the marginal Distributions
    F₁ = cdf.(f₁,rand(f₁,res))
    F₂ = cdf.(f₂,rand(f₂,res))
    data = hcat(F₁, F₂)'

    # Fit the selected copula
    copula_fit = fit(CopType, data)
    # Corresponding joint probability distribution
    joint_pdf = SklarDist(copula_fit,(f₁,f₂))  

    # Evaluation of copula & plotting
    ## Grid of the unit cube [0,1]² for plotting
    x_grid = range(0, stop=1, length=grid_size)
    y_grid = range(0, stop=1, length=grid_size)
    ## Grid on the original scale
    x_grid_OG = range(minimum(X₁), stop=maximum(X₁), length=grid_size)
    y_grid_OG = range(minimum(X₂), stop=maximum(X₂), length=grid_size)
    ## Transform the OG scaled grid to the probability scale using the CDFs
    F₁OG = cdf.(f₁,x_grid_OG)
    F₂OG = cdf.(f₂,y_grid_OG)

    # Evaluate the copula density on the grid
    dens = zeros(grid_size, grid_size)
    dens_OG = zeros(grid_size, grid_size)
    
    for i in 1:grid_size
        for j in 1:grid_size
            # Create a point as a vector of length 2
            point = [x_grid[i], y_grid[j]]  # Point on the unit cube grid
            OG_point = [F₁OG[i], F₂OG[j]]   # Point on the OG scaled grid transformed to probability scale
            # Evaluate the PDF at this point
            dens[i, j] = pdf(copula_fit, point) # Probability density of point on unit cube grid
            dens_OG[i, j] = pdf(copula_fit, OG_point) # Probability density corresponding to OG point
        end
    end

    ## Plot the copula densities
    plt_cop = contour(x_grid, y_grid, dens, fill=true, c=:viridis, title="Copula Density")
    plt_cop_OG = contour(x_grid_OG, y_grid_OG, dens_OG,fill=true, c=:viridis, title="Copula Density - Original Scale")

    return joint_pdf, copula_fit, plt_cop, plt_cop_OG
end

function copula3d_fit_eval(CopType,X₁,X₂,X₃,f₁,f₂,f₃,res,grid_size)
    # X₁, X₂, X₃:     original samples
    # f₁, f₂, f₃:     pdfs of marginal distributions
    # res:        resolution of sampling from pdfs
    # grid_size:  resolution of the grid for the evaluation of the copula
    
    # Calculate the CDFs of the marginal Distributions
    F₁ = cdf.(f₁,rand(f₁,res))
    F₂ = cdf.(f₂,rand(f₂,res))
    F₃ = cdf.(f₃,rand(f₃,res))
    data = hcat(F₁, F₂, F₃)'

    # Fit the selected copula
    copula_fit = fit(CopType, data)
    # Corresponding joint probability distribution
    joint_pdf = SklarDist(copula_fit,(f₁,f₂,f₃))  

    # Evaluation of copula & plotting
    ## Grid of the unit cube [0,1]² for plotting
    x_grid = range(0, stop=1, length=grid_size)
    y_grid = range(0, stop=1, length=grid_size)
    z_grid = range(0, stop=1, length=grid_size)
    ## Grid on the original scale
    x_grid_OG = range(minimum(X₁), stop=maximum(X₁), length=grid_size)
    y_grid_OG = range(minimum(X₂), stop=maximum(X₂), length=grid_size)
    z_grid_OG = range(minimum(X₃), stop=maximum(X₃), length=grid_size)
    ## Transform the OG scaled grid to the probability scale using the CDFs
    F₁OG = cdf.(f₁,x_grid_OG)
    F₂OG = cdf.(f₂,y_grid_OG)
    F₃OG = cdf.(f₃,z_grid_OG)

    # Evaluate the copula density on the grid
    dens = zeros(grid_size, grid_size, grid_size)
    dens_OG = zeros(grid_size, grid_size, grid_size)
    
    for i in 1:grid_size
        for j in 1:grid_size
            for k ∈ 1:grid_size
                # Create a point as a vector of length 2
                point = [x_grid[i], y_grid[j], z_grid[k]]  # Point on the unit cube grid
                OG_point = [F₁OG[i], F₂OG[j], F₃OG[k]]   # Point on the OG scaled grid transformed to probability scale
                # Evaluate the PDF at this point
                dens[i, j, k] = pdf(copula_fit, point) # Probability density of point on unit cube grid
                dens_OG[i, j, k] = pdf(copula_fit, OG_point) # Probability density corresponding to OG point
            end
        end
    end

    ## Plot the copula densities
    plt_cop_xy = contour(x_grid, y_grid, dens[:,:,50], fill=true, c=:viridis, title="Copula Density")
    plt_cop_xz = contour(x_grid, z_grid, dens[:,50,:], fill=true, c=:viridis, title="Copula Density")
    plt_cop_yz = contour(y_grid, z_grid, dens[50,:,:], fill=true, c=:viridis, title="Copula Density")
    plt_cop_OG_xy = contour(x_grid_OG, y_grid_OG, dens_OG[:,:,50],fill=true, c=:viridis, title="Copula Density - Original Scale")
    plt_cop_OG_xz = contour(x_grid_OG, z_grid_OG, dens_OG[:,50,:],fill=true, c=:viridis, title="Copula Density - Original Scale")
    plt_cop_OG_yz = contour(y_grid_OG, z_grid_OG, dens_OG[50,:,:],fill=true, c=:viridis, title="Copula Density - Original Scale")

    return joint_pdf, copula_fit, plt_cop_xy, plt_cop_OG_xy, plt_cop_xz, plt_cop_OG_xz, plt_cop_yz, plt_cop_OG_yz
end

function student_fit(SampVar, initial_params, lower_bounds, upper_bounds)
    # Fit Student's T-distribution
    # initial_params = [μ, σ, ν] : Initial parameter guesses
    # Define bounds for the parameters
    # lower_bounds = [mean(SampVar)-10*std(SampVar); 0.0; 2.0]
    # upper_bounds = [mean(SampVar)+10*std(SampVar); 20.0; 300.0]

    # Optimize the parameters
    results = optimize(negative_log_likelihood, lower_bounds, upper_bounds, initial_params, Fminbox(NelderMead()))

    # Extract the fitted parameters
    fitted_params = Optim.minimizer(results)
    println("t-dist fit parameters: ", fitted_params)
    μ_fit, σ_fit, ν_fit = fitted_params

    fitted_dist = LocationScale(μ_fit, σ_fit, TDist(ν_fit))
    points = range(minimum(SampVar), stop=maximum(SampVar), length=100)
    pdf_values = pdf.(fitted_dist, points)

    plt_hist = histogram(SampVar, normalize=:pdf, lab=:false, palette=[cb[8]; cb[11]; cb[5]])
    plot!(plt_hist, points, pdf_values, label="PDF of fitted t-distribution", xlabel="x", ylabel="Density", title="Fitted Student's t-Distribution PDF")
    
    return plt_hist, fitted_params, fitted_dist
end

# Define the negative log-likelihood function
function negative_log_likelihood(params)
    μ, σ, ν = params
    if σ <= 0 || ν <= 0
        return Inf # Ensure σ and ν are positive
    end
    dist = LocationScale(μ, σ, TDist(ν))
    return -sum(logpdf(dist, x) for x in t)
end

function conditional_sampling(joint_pdf,N,x₁ᶜ,w)
    # For a given joint pdf P(x₁,x₂), sample values of x₂ for given arguments for x₁
    # SampVar: the original dataset
    # joint_pdf: PDF of the joint probability distribution
    # N: number of samples to be drawn

    # Sample from joint probability distribution
    samp_rand = rand(joint_pdf,N)
    x₁_rand = samp_rand[1,:]
    x₂_rand = samp_rand[2,:]

    # Which of the randomly sampled pairs from the joint pdf belong in the specified bin of x₁
    x₁_samp, x₂_samp = [], []
    for i ∈ 1:N
        if x₁_rand[i] < (x₁ᶜ+w/2) && x₁_rand[i] > (x₁ᶜ-w/2)
            push!(x₁_samp,x₁_rand[i])
            push!(x₂_samp,x₂_rand[i])
        end
    end

    x₁_samp, x₂_samp = Float64.(x₁_samp), Float64.(x₂_samp)
    x₁_range = [x₁ᶜ-w/2 x₁ᶜ+w/2]
    # What is the distribution of the x₂ values in the specified x₁ bin
    hist_x2 = histogram(x₂_samp, normalize=:pdf, lab=:false, palette=[cb[8]; cb[11]; cb[5]])
    hist_x1x2 = histogram2d(x₁_samp, x₂_samp, normalize=:pdf, show_empty_bins=true, color=:plasma, xlab = L"\overline{t}", ylab = L"\overline{A}")

    return hist_x2, hist_x1x2, x₁_samp, x₂_samp, x₁_range
end
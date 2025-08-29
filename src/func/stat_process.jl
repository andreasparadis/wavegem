function sample_stats(dataset)
    NoSamp = length(dataset)  # Number of samples in dataset
    stats = zeros(Float64, 1, 6)    # [μ σ med() IQR skew() kurt()]

    stats[1] = mean(dataset)
    stats[2] = std(dataset)
    stats[3] = median(dataset)
    temp_sort = sort(dataset)

    if mod(NoSamp,2) == 0
        Q1 = median(temp_sort[1:Int(NoSamp/2)])
        Q3 = median(temp_sort[Int(NoSamp/2+1):end])
        stats[4] = Q3.-Q1
    else
        Q1 = median(temp_sort[1:Int((NoSamp-1)/2)])
        Q3 = median(temp_sort[Int((NoSamp+3)/2):end])
        stats[4] = Q3.-Q1
    end

    stats[5] = skewness(dataset)
    stats[6] = kurtosis(dataset)

    return stats
end

function sample_stats_multi(datasets)
    NoSamp = size(datasets)[1]  # Number of samples in each dataset
    NoVars = size(datasets)[2]  # Total number of variables
    stats = zeros(Float64, NoVars, 6)    # [μ σ med() IQR skew() kurt()]

    stats[:,1] = mean(datasets, dims=1)' 
    stats[:,2] = std(datasets, dims=1)'
    stats[:,3] = median(datasets, dims=1)'
    temp_sort = sort(datasets, dims=1)

    if mod(NoSamp,2) == 0
        Q1 = median(temp_sort[1:Int(NoSamp/2),:], dims=1)
        Q3 = median(temp_sort[Int(NoSamp/2+1):end,:], dims=1)
        stats[:,4] = Q3.-Q1
    else
        Q1 = median(temp_sort[1:Int((NoSamp-1)/2),:], dims=1)
        Q3 = median(temp_sort[Int((NoSamp+3)/2):end,:], dims=1)
        stats[:,4] = Q3.-Q1
    end
    for i ∈ 1:NoVars
        stats[i,5] = skewness(datasets[:,i])
        stats[i,6] = kurtosis(datasets[:,i])
    end

    return stats
end

function remove_outliers(dataset,cσ)
    # Remove outliers based on the std of the sample
    # cσ: factor of std for threshold of outliers
    L = length(dataset)
    proc_dset = Array{Float64}(undef,0) # Processed dataset (excluding outliers)
    ipd = Array{Int64}(undef,0) # Indices of qualified data points

    stats = sample_stats(dataset)

    for i ∈ 1:L
        if abs(dataset[i]) < stats[3] + cσ*stats[2]
            push!(proc_dset, dataset[i])
            push!(ipd, i)
        end
    end
    
    return proc_dset, ipd
end

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
    # elseif dist_type == Weibull || dist_type == LogNormal
    #     value_range = range(0, maximum(dataset), length=bins)
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
    elseif dist_type == Rayleigh
        param_lab = "$dist_type : σ=$(round(params[1]*1e3)/1e3)"
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
    plt_cop_xy = contour(x_grid, y_grid, dens[:,:,Int(grid_size/2)], fill=true, c=:viridis, title="Copula Density")
    plt_cop_xz = contour(x_grid, z_grid, dens[:,Int(grid_size/2),:], fill=true, c=:viridis, title="Copula Density")
    plt_cop_yz = contour(y_grid, z_grid, dens[Int(grid_size/2),:,:], fill=true, c=:viridis, title="Copula Density")
    plt_cop_OG_xy = contour(x_grid_OG, y_grid_OG, dens_OG[:,:,Int(grid_size/2)],fill=true, c=:viridis, title="Copula Density - Original Scale")
    plt_cop_OG_xz = contour(x_grid_OG, z_grid_OG, dens_OG[:,Int(grid_size/2),:],fill=true, c=:viridis, title="Copula Density - Original Scale")
    plt_cop_OG_yz = contour(y_grid_OG, z_grid_OG, dens_OG[Int(grid_size/2),:,:],fill=true, c=:viridis, title="Copula Density - Original Scale")

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
    hist_x1x2 = histogram2d(x₁_samp, x₂_samp, normalize=:pdf, show_empty_bins=true, color=:plasma)

    return x₁_samp, x₂_samp, x₁_range, hist_x2, hist_x1x2
end

function copula_2D_histogram(joint_pdf,N)
    # For a given joint pdf P(x₁,x₂), sample pairs (x₁,x₂) and generate a 2D histogram from the sampled values 
    # joint_pdf: PDF of the joint probability distribution
    # N: number of samples to be drawn

    # Sample from joint probability distribution
    samp_rand = rand(joint_pdf,N)
    x₁_samp = samp_rand[1,:]
    x₂_samp = samp_rand[2,:]

    # What is the distribution of the x₂ values in the specified x₁ bin
    hist_x1x2 = histogram2d(x₁_samp, x₂_samp, normalize=:pdf, show_empty_bins=true, color=:plasma)

    return hist_x1x2, x₁_samp, x₂_samp
end

function copula_2D_hmap(joint_pdf,N)
    # For a given joint pdf P(x₁,x₂), sample pairs (x₁,x₂) and generate a 2D hsitogram from the sampled values 
    # joint_pdf: PDF of the joint probability distribution
    # N: number of samples to be drawn

    # Sample from joint probability distribution
    samp_rand = rand(joint_pdf,N)
    x₁_samp = samp_rand[1,:]
    x₂_samp = samp_rand[2,:]

    # What is the distribution of the x₂ values in the specified x₁ bin
    hist_x1x2,_ = joint_pdf_hmap(x₁_samp, x₂_samp, 1)

    return hist_x1x2, x₁_samp, x₂_samp
end

function copula_3D_histograms(joint_pdf,N)
    # For a given joint pdf P(x₁,x₂), sample triplets (x₁,x₂) and generate three 2D hsitograms from the sampled values 
    # joint_pdf: PDF of the joint probability distribution
    # N: number of samples to be drawn

    # Sample from joint probability distribution
    samp_rand = rand(joint_pdf,N)
    x₁_samp = samp_rand[1,:]
    x₂_samp = samp_rand[2,:]
    x₃_samp = samp_rand[3,:]

    # What is the distribution of the x₂ values in the specified x₁ bin
    hist_x1x2 = histogram2d(x₁_samp, x₂_samp, normalize=:pdf, show_empty_bins=true, color=:plasma)
    hist_x1x3 = histogram2d(x₁_samp, x₃_samp, normalize=:pdf, show_empty_bins=true, color=:plasma)
    hist_x2x3 = histogram2d(x₂_samp, x₃_samp, normalize=:pdf, show_empty_bins=true, color=:plasma)

    return hist_x1x2, hist_x1x3, hist_x2x3, x₁_samp, x₂_samp, x₃_samp
end

function joint_pdf_hmap(x₁, x₂, fbin)
    # fbin: Bin flag = 0: sqrt, [1:∞]: Freedman-Diaconis rule
    n₁ = length(x₁)
    n₂ = length(x₂)

    x₁ᵘᵇ = maximum(x₁) 
    x₁ˡᵇ = minimum(x₁) 
    x₂ᵘᵇ = maximum(x₂) 
    x₂ˡᵇ = minimum(x₂) 

    if fbin == 0
        # Number of bins
        N₁ = Int(floor(sqrt(n₁))+1)
        N₂ = Int(floor(sqrt(n₂))+1)
    else 
        temp_sort = sort(x₁)
        if mod(n₁,2) == 0
            Q1 = median(temp_sort[1:Int(n₁/2)])
            Q3 = median(temp_sort[Int(n₁/2+1):end])
            IQR = Q3.-Q1
        else
            Q1 = median(temp_sort[1:Int((n₁-1)/2)])
            Q3 = median(temp_sort[Int((n₁+3)/2):end])
            IQR = Q3.-Q1
        end
        wᵇ = 2*IQR*n₁^(-1/3)
        N₁ = Int(floor((x₁ᵘᵇ-x₁ˡᵇ)/wᵇ)) + 1  

        temp_sort = sort(x₂)
        if mod(n₂,2) == 0
            Q1 = median(temp_sort[1:Int(n₂/2)])
            Q3 = median(temp_sort[Int(n₂/2+1):end])
            IQR = Q3.-Q1
        else
            Q1 = median(temp_sort[1:Int((n₂-1)/2)])
            Q3 = median(temp_sort[Int((n₂+3)/2):end])
            IQR = Q3.-Q1
        end
        wᵇ = 2*IQR*n₂^(-1/3)
        N₂ = Int(floor((x₂ᵘᵇ-x₂ˡᵇ)/wᵇ)) + 1 
    end

    x₁ᵇ = range(x₁ˡᵇ, x₁ᵘᵇ, length=N₁)
    x₂ᵇ = range(x₂ˡᵇ, x₂ᵘᵇ, length=N₂)

    hist = fit(Histogram, (x₁, x₂), (x₁ᵇ,x₂ᵇ), closed=:right)
    e₁ = hist.edges[1]
    e₂ = hist.edges[2]
    Wᴴ = hist.weights
    probabs = (Wᴴ / sum(Wᴴ))'
    hmap_x1x2 = contourf(e₁[1:end-1],e₂[1:end-1],probabs*100, color=:plasma, colorbar_title="P[%]",lw=0)

    return hmap_x1x2, probabs, e₁, e₂, Wᴴ
end

function distros_from_jpd(x₁, x₂, dist_type,leg_var)
    # x₁: Dataset of conditioning variable
    # x₂: Dataset of contidional variable
    fbin = 1    # Bin flag = 0: sqrt, [1:∞]: Freedman-Diaconis rule
    Nˢₘᵢₙ = 80 # Minimum number of samples in bin for distribution fitting
    L = length(x₁)  # Number of samples in datasets

    hmap, P₁₂, e₁, e₂, W₁₂ = joint_pdf_hmap(x₁, x₂, fbin)

    Nᵇ₁ = length(e₁)-1 # Number of bins for x₁
    samples = zeros(Int64,Nᵇ₁)  # Number of samples per bin
    [samples[i] = sum(W₁₂[i,:]) for i ∈ 1:Nᵇ₁]

    # x₁ Lower and upper bounds
    iˡᵇ = 1
    temp_iˡᵇ = findall(diff(sign.(samples.-Nˢₘᵢₙ)) .== 2)
    if temp_iˡᵇ != Int64[]
        iˡᵇ = temp_iˡᵇ[1]
    end

    iᵘᵇ = Nᵇ₁
    temp_iᵘᵇ = findall(diff(sign.(samples.-Nˢₘᵢₙ)) .== -2)
    if temp_iᵘᵇ != Int64[]
        iᵘᵇ = temp_iᵘᵇ[1]
    end
    
    # Final bin edges
    x₁ᵇ = e₁[iˡᵇ:iᵘᵇ]
    # Bin widths
    w₁ᵇ = e₁[2] - e₁[1]
    # Bin central values
    x₁ᵇᶜ = x₁ᵇ[1:end-1] .+ w₁ᵇ/2
    Nᴮ = length(x₁ᵇᶜ)

    # Compute the distribution in each x₁ bin based on the x₂ samples that belong to it
    ## Initialisation
    x₂_Fits = zeros(Float64,Nᴮ,2)
    x₂_Modes = zeros(Float64,Nᴮ)
    plt_bin_distros = plot(title="Distribution per bin", palette=:darkrainbow)
    Ø = Array{Float64}(undef,0)

    for i ∈ 1:Nᴮ
        x₂_samp = Ø[:]
        for j ∈ 1:L
            if x₁[j] < x₁ᵇ[i+1] && x₁[j] > x₁ᵇ[i]
                push!(x₂_samp,x₂[j])
            end
        end
        println("Number of samples in bin x₁= $(round(x₁ᵇᶜ[i]*1000)/1000): $(length(x₂_samp))")

        Nbins = 2^9
        pdf_values, fitted_dist, dist_params, _, _, _, x₂_range = dist_fit(x₂_samp, dist_type, Nbins)
        plot!(plt_bin_distros, x₂_range, pdf_values, lw=2, lab=leg_var*"=$(round(x₁ᵇᶜ[i]*1000)/1000)")

        x₂_Fits[i,:] = dist_params[:]
        x₂_Modes[i] = modes(fitted_dist)[1]

        plot!(hmap, [x₁ᵇ[i]; x₁ᵇ[i]], [e₂[1]; e₂[end]], ls=:dash, leg=:false)
        plot!(hmap, [x₁ᵇ[i+1]; x₁ᵇ[i+1]], [e₂[1]; e₂[end]], ls=:dash, leg=:false)
    end

    return x₁ᵇᶜ, x₂_Fits, x₂_Modes, hmap, plt_bin_distros
end

function read_ofast_stats(case_str,case_id,FirstRun,NoRuns)
    # Read from OpenFAST simulation results
    NoEv = zeros(Int8,NoRuns)   # Number of events per run
    TotNoEv, RunCnt = 0, 0
    FairTen = []
    CaseStats = zeros(Float64,NoRuns,6)
    Hₛᴿ = zeros(Float64,NoRuns)
    Tₚᴿ = zeros(Float64,NoRuns)

    for run_id ∈ FirstRun:(FirstRun+NoRuns-1)
        RunCnt += 1
        # Paths
        postOFpath = joinpath(libpath,"SE","$(case_id)","$(run_id)","0","postOFAST")
        fid_tinsts = joinpath(postOFpath,"MaxFair_tinsts_1ptile")
        fid_mfstat = joinpath(postOFpath,case_str,"StatMetrics")
        Decpath = joinpath(libpath,"SE","$(case_id)","$(run_id)","Decomposition")
        fsimpar = joinpath(Decpath,"sim_pars")       # Sea-state parameters file path

        # Parse files
        instances = parse_fxw(fid_tinsts, 0)
        CaseStats[RunCnt,:] = parse_fxw(fid_mfstat, 1) # CaseStats = [μ σ med() rms() skew() kurt()]
        SimPars = parse_fxw(fsimpar,1)
        Hₛᴿ[RunCnt] = SimPars[1]
        Tₚᴿ[RunCnt] = SimPars[2]

        NoEv[RunCnt] = size(instances)[1]
        TotNoEv += NoEv[RunCnt]
        for i ∈ 1:NoEv[RunCnt]
            push!(FairTen, instances[i,2])
        end
    end

    FairTen = Float64.(FairTen);
    H̃ₛᴿ = median(Hₛᴿ) # Median significant wave height from all runs
    T̃ₚᴿ = median(Tₚᴿ) # Median peak wave period from all runs
    MaxNoEv = maximum(NoEv) # Highest number of events in all Runs - For matrix initialisation

    # Candlestick plot for Fairled tensions per run
    plt_CStats = plot(legend=false, title="Fairlead tension response", xlabel="Simulation", ylabel=L"\mu \pm \sigma [N]")
    for i in 1:NoRuns
        plot!(plt_CStats, [i, i], [CaseStats[i,1], CaseStats[i,1]], linewidth=3) # mean
        plot!(plt_CStats, [i, i], [CaseStats[i,1] - CaseStats[i,2], CaseStats[i,1] + CaseStats[i,2]], linewidth=1) # ±std
    end
    plot!(plt_CStats)

    # Effective steepness
    ϵₛₛ = Hₛᴿ ./ Tₚᴿ.^2 * 2π/g

    return FairTen, Hₛᴿ, Tₚᴿ, ϵₛₛ, H̃ₛᴿ, T̃ₚᴿ, NoEv, TotNoEv, plt_CStats
end

function read_ergees_pars(case_str,case_id, FirstRun, NoRuns, TotNoEv, t̅)
    # Read and assign event and Gaussian Elementary Envelopes' (GEEs) parameters

    # Inputs
    Nₜ = length(t̅)
    Nₛ₂ = Int(nextpow(2,Nₜ)/2)
    Ø = Array{Float64}(undef,0)

    ## Initialise variables 
    allΩ, allωₚ, allβ̇, allT₂₋, allΔTₑᵥ, allN, maxAₒ  = [Ø[:] for _ = 1:7]   # For event parameters - Length TotNoEv
    allβ, alltᶜ, alltᶜᵣ, alltᶜₑ, allAₒ, allTₒ, allTₒᵣ, allTₒₑ = [Ø[:] for _ = 1:8] # For GEEs' parameters - Length TotNoGEEs
    allΔtᶜ, allΔtᶜᵣ, allΔtᶜₑ = [Ø[:] for _ = 1:3]    # Time difference between t̅ᶜ of each GEE - Length < TotNoGEEs
    T̅metric = Ø[:]

    ### For means
    allG = zeros(Float64,Nₜ,TotNoEv)    # For event envelopes G(t)
    all_sp2, all_phi = zeros(Float64,Nₛ₂,TotNoEv), zeros(Float64,Nₛ₂,TotNoEv)   # For means of 2nd-
    fsp2 = range(0.0,4.0,Nₛ₂)   # Frequency range
    TotEvID = 0 # Temporary counter for events

    # Initialise plots
    pltG = plot(xlab = L"\overline{t}", ylab = L"\overline{A}", palette=:darkrainbow)
    pltGEE1 = plot3d(xlab = L"\overline{T}", ylab = L"\overline{t}_c", zlab = L"\overline{A}", legend=:false)
    pltGEE2 = plot(xlab = L"\overline{T}", ylab = L"\overline{A}", legend=:false)
    pltGEE3 = plot(xlab = L"\overline{t}_c", ylab = L"\overline{T}", legend=:false)
    pltGEE4 = plot(xlab = L"\overline{t}_c", ylab = L"\overline{A}", legend=:false)
    pltGEE5 = plot(xlab = L"\overline{t}_c", ylab = L"\tilde{\beta}"*" [rad]", legend=:false)
    pltGEE6 = plot(xlab = "Event", ylab = L"\Delta\overline{t}_c", legend=:false)
    pltSP2 = plot(xlab = "f [Hz]", ylab = L"S [m^2 s]", palette=:darkrainbow)
    pltPHI2 = plot(xlab = "f [Hz]", ylab = L"\phi [rad]", palette=:darkrainbow)

    # Read through text files from all runs
    for run_id ∈ FirstRun:(FirstRun+NoRuns-1)
        Decpath = joinpath(libpath,"SE","$(case_id)","$(run_id)","Decomposition")   # Path to Decomposition folder
        for evID ∈ 1:NoEv[run_id-FirstRun+1]
            # Paths
            evdir = joinpath(case_str,"EV$evID")
            Eventpath = joinpath(Decpath,evdir)
            GRrespath = joinpath(Eventpath,"ERGEEs")
            fid_pEWG = joinpath(GRrespath,"EWG_norm_pars")
            fid_pG = joinpath(GRrespath,"G_pars")
            fid_evpars = joinpath(Eventpath,"ev_pars")
            fid_eta2nd = joinpath(Eventpath,"eta_2nd-")
            # fid_eta2nd = joinpath(Eventpath,"event_lin")

            # Read EWGs parameters
            cont = parse_fxw(fid_pEWG, 1)   
            t̅ᶜ, T̅ₒ, A̅ₒ, β̃ = cont[:,1], cont[:,2], cont[:,3], cont[:,4]
            lenT = length(t̅ᶜ)
            cont2 = parse_fxw(fid_pG, 1)     # Read event and GR parameters
            Hₛ, Tₚ, Ω = cont2[1:3]
            ϵₛₛ = Hₛ ./ Tₚ.^2 * 2π/g  # Effective steepness

            # Read parameters of Gaussian Regression of wave event
            cont0 = parse_fxw(fid_evpars, 1)     # Read event parameters
            T₂₋, ΔTₑᵥ = cont0[4], cont0[5]

            ΔT̅ₑᵥ = ΔTₑᵥ/(ϵₛₛ*T₀ₛ)   # Normalisation of ΔTₑᵥ
            push!(allT₂₋,T₂₋)
            push!(allΔTₑᵥ, ΔT̅ₑᵥ)
            # push!(allN, lenT/ΔT̅ₑᵥ)
            push!(allN, lenT/(ΔTₑᵥ/T₀ₛ))

            push!(allΩ, Ω)
            push!(allωₚ, 2π/Tₚ)
            push!(maxAₒ, maximum(A̅ₒ))

            ΣT̅ = 0.0
            for EWG ∈ 1:lenT
                push!(allAₒ,A̅ₒ[EWG])
                push!(allTₒ,T̅ₒ[EWG])
                push!(alltᶜ,t̅ᶜ[EWG])
                push!(alltᶜᵣ,t̅ᶜ[EWG]*Tₚ/T₀ₚ)
                push!(alltᶜₑ,t̅ᶜ[EWG]*Tₚ/nextpow(2,ΔTₑᵥ))
                push!(allβ,β̃[EWG])
                push!(allTₒᵣ,T̅ₒ[EWG]*Tₚ/T₀ₚ)
                push!(allTₒₑ,T̅ₒ[EWG]*Tₚ/nextpow(2,ΔTₑᵥ))
                ΣT̅ += T̅ₒ[EWG]
            end

            # Peak instance intervals
            s_t̅ᶜ = sort(t̅ᶜ)
            for EWG ∈ 1:lenT-1
                push!(allΔtᶜ,s_t̅ᶜ[EWG+1]-s_t̅ᶜ[EWG])
            end

            ## Re-normalise t̅ᶜ to the FOWT dynamics for the Δtᶜ calculation
            # t̅ᶜᵣ = Tₚ*t̅ᶜ./T₀ₛ
            t̅ᶜᵣ = t̅ᶜ.*(Tₚ/T₀ₚ)
            s_t̅ᶜᵣ = sort(t̅ᶜᵣ)
            for EWG ∈ 1:lenT-1
                push!(allΔtᶜᵣ,s_t̅ᶜᵣ[EWG+1]-s_t̅ᶜᵣ[EWG])
            end

            ## Re-normalise t̅ᶜ to the event duration for the Δtᶜ calculation
            t̅ᶜₑ = t̅ᶜ.*(Tₚ/nextpow(2,ΔTₑᵥ))
            s_t̅ᶜₑ = sort(t̅ᶜₑ)
            for EWG ∈ 1:lenT-1
                push!(allΔtᶜₑ,s_t̅ᶜₑ[EWG+1]-s_t̅ᶜₑ[EWG])
            end

            # linfit_β̃ = linear_fit(t̅ᶜ,β̃)
            linfit_β̃ = linear_fit(t̅ᶜᵣ,β̃)
            # linfit_β̃ = linear_fit(t̅ᶜₑ,β̃)
            if !isnan(linfit_β̃[2])
                push!(allβ̇,linfit_β̃[2])
            end

            T̅metric = push!(T̅metric, ΣT̅*Tₚ/ΔTₑᵥ)

            TotEvID += 1

            # Resulting Gaussian approximation of the envelope
            G̅ = gauss_fun(t̅, A̅ₒ, t̅ᶜ, T̅ₒ)
            allG[:,TotEvID] = G̅[:]

            # Read 2nd minus spectrum
            cont_2nd = parse_fxw(fid_eta2nd, 0)     
            time_ev = cont_2nd[:,1]
            eta2nd = cont_2nd[:,2]
            freq_2nd, mag_2nd,phi_2nd,_ = one_side_asp(eta2nd, time_ev)

            # Re-sampling spectrum
            itp_sp2 = interpolate(freq_2nd, mag_2nd, BSplineOrder(4))
            itp_mag_2nd = itp_sp2.(fsp2)
            # Re-sampling phase
            itp_phi = interpolate(freq_2nd, unwrap(phi_2nd), BSplineOrder(4))
            itp_phi_2nd = itp_phi.(fsp2)

            all_sp2[:, TotEvID] = itp_mag_2nd[:]
            all_phi[:, TotEvID] = itp_phi_2nd[:]

            plot!(pltG, t̅, G̅, lab=:false, line=:dot)
            plot!(pltGEE1, [T̅ₒ], [t̅ᶜ], [A̅ₒ], seriestype=:scatter)
            plot!(pltGEE2, [T̅ₒ], [A̅ₒ], seriestype=:scatter, xlim=(0,xlims(pltGEE1)[2]), ylim=(0,zlims(pltGEE1)[2]))
            plot!(pltGEE3, [t̅ᶜ], [T̅ₒ], seriestype=:scatter)
            plot!(pltGEE4, [t̅ᶜ], [A̅ₒ], seriestype=:scatter)
            plot!(pltGEE5, [t̅ᶜ], [β̃], seriestype=:scatter)
            # plot!(pltGEE5, [t̅ᶜ], [unwrap(β̃)], seriestype=:scatter)
            plot!(pltGEE6, TotEvID*ones(Int64,lenT-1), diff(sort(t̅ᶜ)), seriestype=:scatter)
            # plot!(pltGEE6, [TotEvID], [mean(diff(sort(t̅ᶜ)))], seriestype=:scatter)
            # plot!(pltGEE6, [ΔTₑᵥ], [std(diff(sort(t̅ᶜ)))], seriestype=:scatter)
            plot!(pltSP2, fsp2, itp_mag_2nd, lab=:false, line=:dot)
            plot!(pltPHI2, fsp2, itp_phi_2nd, lab=:false, line=:dot)
        end
    end

    Lβ̇ = length(allβ̇)
    if Lβ̇ ≠ TotNoEv
        allβ̇ = [allβ̇; mean(allβ̇)*ones(TotNoEv-Lβ̇)]
    end

    # Critical Wave Events datasets
    CWE_datasets = [allΩ  allβ̇  allΔTₑᵥ allN maxAₒ]
    # Gaussian Elementary Envelopes' datasets
    GEE_datasets = [allAₒ allTₒ allTₒᵣ allTₒₑ alltᶜ alltᶜᵣ alltᶜₑ allβ]
    # Derivative datasets
    Δt_datasets = [allΔtᶜ allΔtᶜᵣ allΔtᶜₑ]
    # Other datasets
    other_datasets = [allωₚ allT₂₋]

    # Mean envelope from all TotNoEv events
    Gmean = mean!(ones(Nₜ,1),allG[:,1:TotEvID])
    plot!(pltG,t̅, Gmean, lw=3, lab="Mean")
    # Mean spectrum of 2nd- components
    sp2mean = mean!(ones(Nₛ₂,1),all_sp2[:,1:TotEvID])
    sp2mean = sp2mean[:,1]
    plot!(pltSP2, fsp2, sp2mean, lw=3, lab="Mean", xlim=(0,0.5))
    # Mean phase angle of 2nd- components
    phi2mean = mean!(ones(Nₛ₂,1),all_phi[:,1:TotEvID])
    plot!(pltPHI2, fsp2, phi2mean, lw=3, lab="Mean")
    # Peak period of mean 2nd- spectrum
    Tᵖ₂₋ = 1/fsp2[findmax(sp2mean)[2]]

    # Reconstruction of 2nd- component via ifft based on the mean spectrum and phase angle
    H2 = Nₛ₂*sp2mean.*exp.(1im*phi2mean)
    # H2 = Nₛ₂*sp2mean.*exp.(1im*2π*rand(Nₛ₂))
    H2 = vcat(H2,0,conj(H2[end:-1:2]))
    η₂₋ = real(ifft(H2))
    η₂₋ = [η₂₋[Int(Nₜ/2):end]; η₂₋[1:Int(Nₜ/2-1)]]

    plt_eta2 = plot(t̅, η₂₋, xlab = L"\overline{t}", ylab = L"\eta_{2nd} ~[m]", lw=2)

    ## Scatter plots of datasets
    pltΩ = plot([allΩ], seriestype=:scatter, legend=:false ,xlab="Event", ylab=L"\Omega=\frac{\tilde{\omega}}{\omega_p}~[-]")
    pltT2 = plot([allT₂₋], seriestype=:scatter, legend=:false ,xlab="Event", ylab=L"T_{2^-}~[s]")
    pltΔT = plot([allΔTₑᵥ], seriestype=:scatter, legend=:false ,xlab="Event", ylab=L"\overline{\Delta T}_{ev}")
    pltN = plot([allN], seriestype=:scatter, legend=:false ,xlab="Event", ylab=L"\frac{No WGs}{\overline{\Delta T}_{ev}}")

    return CWE_datasets, GEE_datasets, Δt_datasets, other_datasets, Tᵖ₂₋, T̅metric, pltGEE1, pltGEE2, pltGEE3, pltGEE4, pltGEE5, pltGEE6, pltG, pltSP2, pltPHI2, plt_eta2
end

function tc_Atc_dist()
    # Just keeping this script here for future reference.
    # Trying the tᶜ-A*tᶜ distribution
    ## Remove tᶜ=0 values to avoid NaN results and also remove outliers
    alltᶜ_nz = Array{Float64}(undef,0)
    iz = Array{Int64}(undef,0)

    for i ∈ 1:length(alltᶜ)
        if alltᶜ[i] ≠ 0 && abs(alltᶜ[i]) < 3*std(alltᶜ)
            push!(alltᶜ_nz, alltᶜ[i])
            push!(iz, i)
        end
    end

    Aotc = allAₒ[iz]./alltᶜ_nz

    Aotc_bins = range(minimum(Aotc), maximum(Aotc), length=Int(round(sqrt(length(Aotc)))))
    hist_Aotc = fit(Histogram, Aotc, Aotc_bins, closed=:right)
    p_Aotc = (hist_Aotc.weights / sum(hist_Aotc.weights))
    plot(Aotc_bins[1:end-1], p_Aotc)

    μᴬᵗ = mean(abs.(Aotc)); σᴬᵗ = std(abs.(Aotc))

    Aotc_fin, ifin = [], []
    for i ∈ 1:length(Aotc)
        if abs(Aotc[i]) < (Aotc_bins[2]-Aotc_bins[1])/2
            push!(Aotc_fin, Aotc[i])
            push!(ifin, i)
        end
    end

    Aotc_fin = Float64.(Aotc_fin)
    alltᶜ_fin = alltᶜ_nz[ifin]

    plot(alltᶜ_fin,Aotc_fin, seriestype=:scatter)#, ylim=(-1,1))
    plot!(sort(alltᶜ_fin), π/2*exp.(-abs.(sort(alltᶜ_fin))).*coth.(sort(alltᶜ_fin)), lw=2)
    # plot!(sort(alltᶜ_fin), exp(1) * sign.(sort(alltᶜ_fin)) .* exp.(-abs.(sort(alltᶜ_fin))), lw=2)
    # plot!(sort(alltᶜ_fin), sign.(sort(alltᶜ_fin)) .* exp.(-sqrt.(abs.(sort(alltᶜ_fin)))), lw=2)

    plot(alltᶜ,allAₒ, seriestype=:scatter)
    plot!(alltᶜ_fin,Aotc_fin.*alltᶜ_fin, seriestype=:scatter)

    histogram2d(alltᶜ_fin,Aotc_fin, normalize=:pdf, show_empty_bins=true, color=:plasma, bins=:sqrt)
    histogram(Aotc_fin, normalize=:pdf, show_empty_bins=true, color=:plasma, bins=:sqrt)

    dataset = Aotc_fin;   dist_type = Cauchy;  Nbins = 2^9
    pdf_values, Aotc_Fit_fin, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
    plot!(hst_plt,  title="PDF - P(A̅ₒ/t̅ᶜ)",  xlab=L"\frac{\overline{A}_o}{\overline{t}_c}", ylab=L"P(\frac{\overline{A}_o}{\overline{t}_c})")
    FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])

    dataset = alltᶜ_fin;   dist_type = Normal;  Nbins = 2^9
    pdf_values, tᶜ_Fit_fin, dist_params, samp_params, hst_plt, QQplt = dist_fit(dataset, dist_type, Nbins)
    plot!(hst_plt,  title="PDF - P(t̅ᶜ)", xlab=L"\overline{t}^c", ylab=L"P(\overline{t}^c)")
    FitPlot = plot(hst_plt, QQplt, layout = @layout[a; b])

    JPD_tᶜAtc, Cop_tᶜAtc, plt_cop_tᶜAtc, plt_cop_tᶜAtc_OG = copula2d_fit_eval(GaussianCopula,alltᶜ_fin,Aotc_fin,tᶜ_Fit_fin,Aotc_Fit_fin,2^10,100)
    ### Add labels to the copula density plots
    plot!(plt_cop_tᶜAtc_OG, xlabel=L"\overline{t}^c", ylabel=L"\frac{\overline{A}_o}{\overline{t}_c}")
end

function poly_from_maxima(x_coords::Vector{Float64}, y_coords::Vector{Float64}, zero_crossings::Vector{Float64})
    n = length(x_coords)
    @assert length(y_coords) == n "The number of x-coordinates and y-coordinates must be the same."

    degree = 2n + 1

    # Number of conditions: n for function values, n for first derivatives, n for second derivatives, and 2 for zero crossings
    num_conditions = 3n + 2

    # Initialize the matrix A and vector B
    A = zeros(Float64, num_conditions, degree + 1)
    B = zeros(Float64, num_conditions)

    # Fill the matrix A and vector B
    for i in 1:n
        x_i = x_coords[i]
        y_i = y_coords[i]

        # Function value conditions
        for j in 0:degree
            A[i, j + 1] = x_i^j
        end
        B[i] = y_i

        # Derivative conditions
        for j in 1:degree
            A[n + i, j + 1] = j * x_i^(j - 1)
        end

        # Second derivative conditions (for maxima)
        for j in 2:degree
            A[2n + i, j + 1] = j * (j - 1) * x_i^(j - 2)
        end
        B[2n + i] = -1  # Ensuring the second derivative is negative
    end

    # Fill the matrix A and vector B for zero-crossing conditions
    for (i, x_z) in enumerate(zero_crossings)
        # Zero-crossing conditions: f(x_z) = 0
        for j in 0:degree
            A[3n + i, j + 1] = x_z^j
        end
        B[3n + i] = 0
    end

    # Solve the system of equations
    coefficients = A \ B

    return coefficients
end

function quadratic_with_minimum(P1,P2)
    # We need a quadratic polynomial: ax^2 + bx + c
    # Conditions:
    # 1. f(x1) = y1
    # 2. f(x2) = y2

    x1, y1, x2, y2 = P1[1], P1[2], P2[1], P2[2]
    # Construct the system of equations
    # We need a third condition to determine the three coefficients uniquely.
    # We assume the vertex (minimum) is midway between x1 and x2 for simplicity.
    x_min = (x1 + x2) / 2

    # The derivative at x_min should be zero: f'(x_min) = 0
    # f'(x) = 2ax + b => 2a*x_min + b = 0

    # Construct the matrix A and vector B for the system of equations
    A = [
        x1^2 x1 1
        x2^2 x2 1
        2x_min 1 0
    ]

    B = [y1, y2, 0]

    # Solve the system of equations
    coefficients = A \ B
    a, b, c = coefficients

    # Define the polynomial
    f(x) = a*x^2 + b*x + c

    # Find the minimum point
    x_min_calculated = -b / (2a)
    y_min = f(x_min_calculated)

    return coefficients, x_min_calculated, y_min
end

function fit_quadratic_polynomial(X, Y)
    x1, x2, x3 = X
    y1, y2, y3 = Y

    # We need a quadratic polynomial: ax^2 + bx + c
    # Conditions:
    # 1. f(x1) = 0
    # 2. f(x2) = y2
    # 3. f(x3) = 0
    # 4. f'(x2) = 0 (since x2 is a maximum)

    # Construct the matrix A and vector B for the system of equations
    A = [
        x1^2 x1 1
        x2^2 x2 1
        x3^2 x3 1
        2x2 1 0
    ]

    B = [0, y2, 0, 0]

    # Solve the system of equations
    coefficients = A \ B
    a, b, c = coefficients

    return coefficients
end
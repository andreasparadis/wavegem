using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles
using FFTW, DSP, Statistics, BSplineKit
import ColorSchemes.darkrainbow

gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))
cb = darkrainbow

# Include necessary scripts for functions
include("text_process.jl"); include("signal_processing.jl")
include("peak_detect.jl");  include("num_als.jl")
include("rigid_dyn.jl");    include("hydrodyn.jl")
include("utilities.jl")

##############################################################################
# Flags
facc = true
fcam = true
fprb = true

# Paths for accelerometer, camera and probe measurements
fdate = "2024-10-07"
dname = "EV1_FCSD_05"
fname = "data_0.txt"

parent = joinpath(pwd(),"library","UCL","Model")
accpath = joinpath(parent,"acc",fdate,dname,fname)
campath = joinpath(parent,"cam",fdate)
prbpath = joinpath(parent,"probes",fdate,dname*".txt")
# prbpath = joinpath(parent,"probes",fdate,"EV3_SEWG_00.txt")
accfigs = joinpath(parent,"acc",fdate,dname)
camfigs = joinpath(parent,"cam",fdate,dname)
prbfigs = joinpath(parent,"probes",fdate)

# General
## FOWT eigenfrequencies (theoretical)
λ = 100 # Scale factor (L_full/L_model)
T₀ₛ, T₀ₕ, T₀ₚ = 110.0, 17, 32.5
T₀ₛ = T₀ₛ./ sqrt(λ)
T₀ₕ = T₀ₕ./ sqrt(λ)
T₀ₚ = T₀ₚ./ sqrt(λ)

# Signal processing
fᶜ = 4      # [Hz] Cut-off frequency for low-pass filtering
fₛᵃ = 100   # [Hz] Accelerometer sampling rate
fₛᶜᵃᵐ = 120 # [Hz] Cameras' sampling rate
Nₜₐₚ = 2^8  # Number of taps

###########################################################################
# Read measurements from Bluetooth accelerometer
open(accpath, "r")
content = readdlm(accpath, skipstart=1, '\n')

NoRows = length(content)
NoCols = length(split(content[1])) - 4

parsed_cont = zeros(Float64,NoRows, NoCols)
Time = zeros(Float64,NoRows)
ChipTime = zeros(Float64,NoRows)

# Time: Cumulative calculation from Col1 (HH:MM:SS)
for i ∈ 1:NoRows
    tmp = split(content[i])
    tmp2 = split(tmp[1],":")
    for j ∈ 1:3
        tmp3 = parse(Float64, tmp2[j])
        Time[i] = Time[i] + tmp3*60^(3-j)
    end
end
Time = Time .- Time[1] # Shift to 0

# ChipTime: Cumulative calculation from Col4 (HH:MM:SS)
for i ∈ 1:NoRows
    tmp = split(content[i])
    tmp2 = split(tmp[4],":")
    for j ∈ 1:3
        tmp3 = parse(Float64, tmp2[j])
        ChipTime[i] = ChipTime[i] + tmp3*60^(3-j)
    end
end
ChipTime = ChipTime .- Time[1] # Shift to 0

# Parse rest of data (Col5 to end)
for i ∈ 1:NoRows
    for j ∈ 1:NoCols
        tmp = split(content[i])
        parsed_cont[i,j] = parse(Float64, tmp[j+4])
    end
end

# Low-pass filtering of accelerometer measurements
parsed_cont = low_pass_filter(parsed_cont,fᶜ,fₛᵃ,Nₜₐₚ)

# Isolate data of interest and assign to struct
head_acc = ("Acceleration X","Acceleration Y","Acceleration Z",
"Roll rate","Pitch rate","Yaw rate",
"Roll","Pitch","Yaw")

𝚶⃗ₙₐ = zeros(Float64,NoRows)
Bᵃ = Acc(𝚶⃗ₙₐ[:],𝚶⃗ₙₐ[:],𝚶⃗ₙₐ[:],𝚶⃗ₙₐ[:],𝚶⃗ₙₐ[:],𝚶⃗ₙₐ[:],𝚶⃗ₙₐ[:],𝚶⃗ₙₐ[:],𝚶⃗ₙₐ[:])
# Originally: positive x-axis towards wavemaker, positive z-axis upwards 
# Conform to camera's CCS: Rotation π around z-axis => Change signs of x,y, ϕ, θ
## Accelerations originally in [g] => [m/s²]
Bᵃ.u̇ = parsed_cont[:,1]*9.81
Bᵃ.v̇ = parsed_cont[:,2]*9.81
Bᵃ.ẇ = parsed_cont[:,3]*9.81
## Angular velocities originally in [°/s] => [rad/s]
Bᵃ.ωˣ = -parsed_cont[:,4]*π/180
Bᵃ.ωʸ = -parsed_cont[:,5]*π/180
Bᵃ.ωᶻ = parsed_cont[:,6]*π/180
## Angles originally in [°] => [rad]
Bᵃ.ϕ = -parsed_cont[:,7]*π/180
Bᵃ.θ = -parsed_cont[:,8]*π/180
Bᵃ.ψ = parsed_cont[:,9]*π/180

###########################################################################
# Read measurements from cameras
if fcam
    # PRIOR TO 2024-10-07: camera CS is rotated -π/2 relative to the accelerometer's CS
    # => x = sway, y = surge, z = heave
    # 2024-10-07 onwards: camera CS is rotated by π relative to the accelerometer's CS
    # => x = surge, y = sway, z = heave
    # As the script stands, use only for measurements from 2024-10-07 onwards

    # X,Y,Z: Marker positions in the global CS of the cameras
    fid = joinpath(campath,"$dname-X.txt")
    cont = parse_fxw(fid, 0)
    t, Xᶜᵃᵐ = cont[:,1], cont[:,2:end]
    fid = joinpath(campath,"$dname-Y.txt")
    cont = parse_fxw(fid, 0)
    Yᶜᵃᵐ = cont[:,2:end]
    fid = joinpath(campath,"$dname-Z.txt")
    cont = parse_fxw(fid, 0)
    Zᶜᵃᵐ = cont[:,2:end]

    # Define the Inertial (I) reference frame with z=0 at the water level
    # Used markers 7,8 on the tower axis
    xₜ = (Xᶜᵃᵐ[1,7]+Xᶜᵃᵐ[1,8])/2
    yₜ = (Yᶜᵃᵐ[1,7]+Yᶜᵃᵐ[1,8])/2
    zₜ = (Zᶜᵃᵐ[1,7]+Zᶜᵃᵐ[1,8])/2

    ϕ₀ = atan.((Zᶜᵃᵐ[1,1]-Zᶜᵃᵐ[1,4]) / (Yᶜᵃᵐ[1,1]-Yᶜᵃᵐ[1,4]))
    θ₀ = atan.((Zᶜᵃᵐ[1,1]-Zᶜᵃᵐ[1,4]) / (Xᶜᵃᵐ[1,1]-Xᶜᵃᵐ[1,4]))
    # The vertical distance of markers 7,8 from still water was measured at 165 mm
    rₜ = [xₜ;yₜ;zₜ] .+ 165*[sin(θ₀);sin(ϕ₀);-1]

    # Translate to the Inertial (I) reference frame
    X = Xᶜᵃᵐ .- rₜ[1]
    Y = Yᶜᵃᵐ .- rₜ[2]
    Z = Zᶜᵃᵐ .- rₜ[3]

    # Low-pass filtering of measurements & transformation from [mm] to [m]
    X = low_pass_filter(X.*1e-3,fᶜ,fₛᶜᵃᵐ,Nₜₐₚ)
    Y = low_pass_filter(Y.*1e-3,fᶜ,fₛᶜᵃᵐ,Nₜₐₚ)
    Z = low_pass_filter(Z.*1e-3,fᶜ,fₛᶜᵃᵐ,Nₜₐₚ)

    # The first column in the .txt files is the time vector
    tₑ = t[end]
    Nₜ = length(t)
    dt = tₑ/(Nₜ-1)
    fₛ = 1/dt
    NoMark = size(X)[2]

    # A struct with the time dependent position vectors of all markers
    MRK = Markers(pos(X[:,1],Y[:,1],Z[:,1]), pos(X[:,2],Y[:,2],Z[:,2]), 
    pos(X[:,3],Y[:,3],Z[:,3]), pos(X[:,4],Y[:,4],Z[:,4]), pos(X[:,5],Y[:,5],Z[:,5]), 
    pos(X[:,6],Y[:,6],Z[:,6]), pos(X[:,7],Y[:,7],Z[:,7]), pos(X[:,8],Y[:,8],Z[:,8]),
    pos(X[:,9],Y[:,9],Z[:,9]), pos(X[:,10],Y[:,10],Z[:,10]))

    # Initialization
    𝚶⃗ₙ = zeros(Float64,Nₜ)
    𝚶⃗ₙₓₘ = zeros(Float64,Nₜ,NoMark)

    ## Body reference frame struct
    B = kin(𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],
                𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:])
    ## Inertial reference frame strcut
    I = kin(𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],
                𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:],𝚶⃗ₙₓₘ[:,:])
    # Center of Gravity kinematics
    CG = point_kin(𝚶⃗ₙ[:], 𝚶⃗ₙ[:], 𝚶⃗ₙ[:], 𝚶⃗ₙ[:], 𝚶⃗ₙ[:], 𝚶⃗ₙ[:], 𝚶⃗ₙ[:], 𝚶⃗ₙ[:], 𝚶⃗ₙ[:], 𝚶⃗ₙ[:], 𝚶⃗ₙ[:], 𝚶⃗ₙ[:],
                𝚶⃗ₙ[:], 𝚶⃗ₙ[:], 𝚶⃗ₙ[:], 𝚶⃗ₙ[:], 𝚶⃗ₙ[:], 𝚶⃗ₙ[:])
    x,y,z = [𝚶⃗ₙₓₘ[:,:] for _ = 1:3]

    # Kinematics from marker trajectories
    for i ∈ 1:NoMark
        # Displacements (Motion of markers relative to initial position in the inertial frame)
        x[:,i] = X[:,i].-X[end,i]
        y[:,i] = Y[:,i].-Y[end,i]
        z[:,i] = Z[:,i].-Z[end,i]

        # Inertial (I) reference frame kinematics
        I.x[:,i], I.y[:,i], I.z[:,i] = X[:,i], Y[:,i], Z[:,i] # Marker trajectory
        ## Call function for calculation of kinematics
        Rᴵ, ϕᴵ, θᴵ, ψᴵ, υᴵ, aᴵ, ωᴵ, ω̇ᴵ = kinematics_FOWT(I.x[:,i], I.y[:,i], I.z[:,i], t)
        ## Assign results to reference frame struct
        I.ϕ[:,i], I.θ[:,i], I.ψ[:,i] = unwrap(ϕᴵ), unwrap(θᴵ), unwrap(ψᴵ)
        I.u[:,i], I.v[:,i], I.w[:,i] = υᴵ[:,1], υᴵ[:,2], υᴵ[:,3]
        I.u̇[:,i], I.v̇[:,i], I.ẇ[:,i] = aᴵ[:,1], aᴵ[:,2], aᴵ[:,3]
        I.ωˣ[:,i], I.ωʸ[:,i], I.ωᶻ[:,i] = ωᴵ[:,1], ωᴵ[:,2], ωᴵ[:,3]
        I.ω̇ˣ[:,i], I.ω̇ʸ[:,i], I.ω̇ᶻ[:,i] = ω̇ᴵ[:,1], ω̇ᴵ[:,2], ω̇ᴵ[:,3]

        # Body-fixed (B) refernce frame kinematics
        ## Using one of the markers as a reference point since CG is unknown up to here
        RM = 4  # Reference marker
        ### Position of each marker relative to the reference point
        X₁ᵢ = X[:,i] .- X[:,RM]
        Y₁ᵢ = Y[:,i] .- Y[:,RM]
        Z₁ᵢ = Z[:,i] .- Z[:,RM]
        ## Calculation of kinematics
        Rᴮ, ϕᴮ, θᴮ, ψᴮ, υᴮ, aᴮ, ωᴮ, ω̇ᴮ = kinematics_FOWT(X₁ᵢ, Y₁ᵢ, Z₁ᵢ, t)
        ## Assign results to reference frame struct
        B.u[:,i], B.v[:,i], B.w[:,i] = υᴮ[:,1], υᴮ[:,2], υᴮ[:,3]
        B.u̇[:,i], B.v̇[:,i], B.ẇ[:,i] = aᴮ[:,1], aᴮ[:,2], aᴮ[:,3]
    end

    # For the calculation of body angles and angular velocities, I choose
    # the reference and relative markers based on the corresponding plane
    # for better agreement with the accelerometer
    ## B.ϕ
    X₁ᵢ, Y₁ᵢ, Z₁ᵢ = X[:,8] .- X[:,7], Y[:,8] .- Y[:,7], Z[:,8] .- Z[:,7]
    _, ϕᴮ, _, _, _, _, ωᴮ, ω̇ᴮ = kinematics_FOWT(X₁ᵢ, Y₁ᵢ, Z₁ᵢ, t)
    for i ∈ 1:NoMark
        B.ϕ[:,i] = unwrap(ϕᴮ .- ϕᴮ[end])
        B.ωˣ[:,i] = ωᴮ[:,1]
        B.ω̇ˣ[:,i] = ω̇ᴮ[:,1] 
    end

    ## B.θ
    X₁ᵢ, Y₁ᵢ, Z₁ᵢ = X[:,3] .- X[:,1], Y[:,3] .- Y[:,1], Z[:,3] .- Z[:,1]
    _, _, θᴮ, _, _, _, ωᴮ, ω̇ᴮ = kinematics_FOWT(X₁ᵢ, Y₁ᵢ, Z₁ᵢ, t)
    for i ∈ 1:NoMark
        B.θ[:,i] = unwrap(θᴮ .- θᴮ[end])
        B.ωʸ[:,i] = ωᴮ[:,2]
        B.ω̇ʸ[:,i] = ω̇ᴮ[:,2]
    end

    ## B.ψ
    X₁ᵢ, Y₁ᵢ, Z₁ᵢ = X[:,3] .- X[:,1], Y[:,3] .- Y[:,1], Z[:,3] .- Z[:,1]
    _, _, _, ψᴮ, _, _, ωᴮ, ω̇ᴮ = kinematics_FOWT(X₁ᵢ, Y₁ᵢ, Z₁ᵢ, t)
    for i ∈ 1:NoMark
        B.ψ[:,i] = unwrap(ψᴮ .- ψᴮ[end])
        B.ωᶻ[:,i] = ωᴮ[:,3]
        B.ω̇ᶻ[:,i] = ω̇ᴮ[:,3] 
    end

    # Assign the final angular quantities resulting from the camera measurements
    ϕ, θ, ψ, ωˣ, ωʸ, ωᶻ, ω̇ˣ, ω̇ʸ, ω̇ᶻ = B.ϕ[:,1], B.θ[:,1], B.ψ[:,1], B.ωˣ[:,1], B.ωʸ[:,1], B.ωᶻ[:,1], B.ω̇ˣ[:,1], B.ω̇ʸ[:,1], B.ω̇ᶻ[:,1]

    # Initial angles of markers in the Inertial frame
    Φ₀ᴵ = I.ϕ[1,:]  
    Θ₀ᴵ = I.θ[1,:]
    Ψ₀ᴵ = I.ψ[1,:]

    # Translation Velocity (For each component, each column should be the same)
    u₀ = mean!(ones(Nₜ), I.u .- B.u)
    v₀ = mean!(ones(Nₜ), I.v .- B.v)
    w₀ = mean!(ones(Nₜ), I.w .- B.w)
    # Translation Acceleration 
    u̇₀ = mean!(ones(Nₜ), I.u̇ .- B.u̇)
    v̇₀ = mean!(ones(Nₜ), I.v̇ .- B.v̇)
    ẇ₀ = mean!(ones(Nₜ), I.ẇ .- B.ẇ)

    # Translational motion
    x₀,y₀,z₀ = [𝚶⃗ₙ[:] for _ = 1:3]
    for i ∈ 2:Nₜ
        x₀[i] = (u₀[i]+u₀[i-1])/2 * (t[i]-t[i-1]) + x₀[i-1]
        y₀[i] = (v₀[i]+v₀[i-1])/2 * (t[i]-t[i-1]) + y₀[i-1]
        z₀[i] = (w₀[i]+w₀[i-1])/2 * (t[i]-t[i-1]) + z₀[i-1]
    end

    # Center of gravity 
    zG0 = -163.31*1e-3 # Calculated from function (in rigid_dyn.jl now)
    xG0 = zG0*sin(θ₀)
    yG0 = zG0*sin(ϕ₀)
    rᵢᶜᵍ = [xG0; yG0; zG0]

    CG.x = xG0 .+ x₀
    CG.y = yG0 .+ y₀
    CG.z = zG0 .+ z₀
    CG.u, CG.v, CG.w = u₀[:], v₀[:], w₀[:]
    CG.u̇, CG.v̇, CG.ẇ = u̇₀[:], v̇₀[:], ẇ₀[:]

    # CG acceleration - Use for comparison with accelerometer
    u̇₀,_ = derivatives(u₀,t)
    v̇₀,_ = derivatives(v₀,t)
    ẇ₀,_ = derivatives(w₀,t)

    plot3d(CG.x, CG.y, CG.z, label="Trajectory", linewidth=2)
    scatter3d!([xG0], [yG0], [zG0], label="Initial position", color=:red, markersize=2)
    xlabel!("X")
    ylabel!("Y")
    zlabel!("Z")
    title!("CG - 3D Trajectory")

    # @gif for i ∈ 1:10
    #     plot([CG.x[i]], [CG.z[i]])
    # end

    # Initial position of markers with respect to CG (i.e. in the Body-fixed system)
    rᴮ₀ = [(X[1,:].-xG0) (Y[1,:].-yG0) (Z[1,:].-zG0)]'
    Rᴮ = sqrt.((X[1,:].-xG0).^2 .+ (Y[1,:].-yG0).^2 .+ (Z[1,:].-zG0).^2)

    # Stability (Metacenter - zM must be above zG)
    zM0, zB = float_stabil(0.250,0.100,0.175,0.030)
    xM0 = zM0*sin(θ₀)
    yM0 = zM0*sin(ϕ₀)

    # Sychronize accelerometer & cameras
    ## Find the first peak in corresponding signals and align based on the timing of that peak 
    ### Accelerometer pitch signal
    tₐ = range(0,Time[end],length(Time))    # Accelerometer time
    MinPeakVal, MinPeakDist = std(Bᵃ.θ), 0  # Using standard deviation of signal as threshold value
    PeakVal, PeakPos, PeakId, Der1, Der2 = peaks_extend(Bᵃ.θ,tₐ, MinPeakVal, MinPeakDist)
    τₐ = PeakPos[1] # Time of peak on accelerometer signal

    plt_accpeaks = plot(tₐ, Bᵃ.θ, xlab = "Time [sec]", ylab = "Pitch [deg]",  lab = "Original signal")
    plot!(PeakPos, PeakVal, seriestype=:scatter, ms=2, mc=:red, lab = "Peaks")
    display(plt_accpeaks)
    
    ### Camera pitch signal
    MinPeakVal, MinPeakDist = std(B.θ[:,1]), 0
    PeakVal, PeakPos, PeakId, Der1, Der2 = peaks_extend(B.θ[:,1],t, MinPeakVal, MinPeakDist)
    τᶜ = PeakPos[1] # Time of peak on camera signal

    plt_campeaks = plot(t, B.θ[:,1], xlab = "Time [sec]", ylab = "Pitch [rad/s]",  lab = "Original signal")
    plot!(PeakPos, PeakVal, seriestype=:scatter, ms=2, mc=:red, lab = "Peaks")
    display(plt_campeaks)
    # !!! CHECK THE PLOTS SO THAT PEAK SELECTION WAS DONE CORRECTLY !!!

    Δτ = τₐ - τᶜ # Accelerometer delay
    tₐ = range(-Δτ,Time[end]-Δτ,length(Time))   # Redefine accelerometer time
    
    # Resample the accelerometer signal using the camera's sampling frequency
    itp = interpolate(tₐ, -parsed_cont[:,4]*π/180, BSplineOrder(4))
    ωˣₐ = itp.(t)
    itp = interpolate(tₐ, -parsed_cont[:,5]*π/180, BSplineOrder(4))
    ωʸₐ = itp.(t)
    itp = interpolate(tₐ, parsed_cont[:,6]*π/180, BSplineOrder(4))
    ωᶻₐ = itp.(t)
    itp = interpolate(tₐ, -parsed_cont[:,7]*π/180, BSplineOrder(4))
    ϕₐ = itp.(t)
    itp = interpolate(tₐ, -parsed_cont[:,8]*π/180, BSplineOrder(4))
    θₐ = itp.(t)
    itp = interpolate(tₐ, parsed_cont[:,9]*π/180, BSplineOrder(4))
    ψₐ = itp.(t)

    # Rotational motion of given point with respect to CG
    rᴵₚ = [0.124; 0.124; -0.235] # Initial position vector in inertial reference frame
    rᴮₚ = rᴵₚ - rᵢᶜᵍ # Initial position vector in Body frame
    Φ₀ₚ = atan.(rᴮₚ[3],rᴮₚ[2]) 
    Θ₀ₚ = atan.(rᴮₚ[1],rᴮₚ[3]) 
    Ψ₀ₚ = atan.(rᴮₚ[2],rᴮₚ[1])

    xᵩ = sqrt.(rᴮₚ[1]^2 + rᴮₚ[3]^2)*sin.(θ.+Θ₀ₚ)
    yᵩ = sqrt.(rᴮₚ[2]^2 + rᴮₚ[3]^2)*cos.(ϕ.+Φ₀ₚ)
    zᵩ = sqrt.(rᴮₚ[1]^2 + rᴮₚ[3]^2)*cos.(θ.+Θ₀ₚ)

    # Fairlead tension
    kᵣ = 6.3 # [N/m] Mooring line stiffness
    Fᵢ = kᵣ * (0.6-0.52)
    dr = sqrt.((x₀ .+ xᵩ).^2 .+ (y₀ .+ yᵩ).^2 .+ (z₀ .+ zᵩ).^2)
    Fᶠ = kᵣ*dr .+ Fᵢ

    plt_fairten = plot(t, Fᶠ, xlab="t [s]", ylab="F[N]", lab="Line 1", lw =2)

    # Spectral analysis: 
    ## Surge
    FRx, MAGx, ϕx, _, _, _, Hc = one_side_asp(x₀, t)
    MXₚ = findmax(MAGx)[1]; fˣ = FRx[findmax(MAGx)[2]]
    ## Heave
    FRz, MAGz, ϕz, _, _, _, Hc = one_side_asp(z₀, t)
    MZₚ = findmax(MAGz)[1]; fᶻ = FRz[findmax(MAGz)[2]]
    ## Pitch
    FRθ, MAGθ, ϕθ, _, _, _, Hc = one_side_asp(θ*180/π, t)
    MΘₚ = findmax(MAGθ)[1]; fθ = FRθ[findmax(MAGθ)[2]]

    R = 100*1e-3
    xᵣ = R*cos.(range(0,2π,360)) .- 165*1e-3*sin(θ₀)#.+ r₀[1]
    yᵣ = R*sin.(range(0,2π,360)) .- 165*1e-3*sin(ϕ₀)#.+ r₀[2]
    zᵣ = zeros(Float64,360)
end

###########################################################################
# Read measurements from wave probes
fₚ = 1 # [Hz] Wave peak frequency

if fprb
    prb_cont = parse_fxw(prbpath, 1)

    tₚᵣ = prb_cont[:,1]
    tₚᵣ = tₚᵣ .- tₚᵣ[1]

    Vres = 0.02 # [m/Volt]
    Zₚᵣ = Vres * prb_cont[:,6]
end
###########################################################################
# PLOTS
## Accelerometer
for i ∈ 1:9
    plti = plot(xlab = "t [s]", title = "$(head_acc[i])", legend=:topleft, palette=[cb[11]])
    if i==1 || i==2 || i==3
        plot!(Time, parsed_cont[:,i]*9.81, lw=2, lab = fname, ylab=L"[m/s^2]")
    elseif i==4 || i==5 || i==6
        plot!(Time, parsed_cont[:,i]*π/180, lw=2, lab = fname, ylab=L"[rad/s]")
    else
        plot!(Time, parsed_cont[:,i], lw=2, lab = fname, ylab="[°]")
    end
    display(plti)
    figid = joinpath(accfigs,"$(head_acc[i])")
    savefig(plti,figid*".svg")
    savefig(plti,figid*".png")
end

# Spectral analysis: Acc X (Surge)
FR, MAG, phi, df, T, Nfft, H = one_side_asp(parsed_cont[:,1], Time)
MAGₚ = findmax(MAG)[1]
# fₚ = FR[findmax(MAG)[2]]

# VMAG = zeros(Float64, length(FR))
# DMAG = zeros(Float64, length(FR))
# VMAG[2:end] = abs.(MAG[2:end] ./ (1im * 2π*FR[2:end]))
# DMAG[2:end] = abs.(-MAG[2:end] ./ (2π*FR[2:end]).^2)
# H2 = Nfft*DMAG.*exp.(1im*phi)
# H2 = vcat(H2,0,conj(H2[end:-1:2]))
# surge = ifft(DMAG)

plt_SPX = plot(xlab=L"f~[Hz]", ylab=L"Magnitude~[m/s^2]", palette=[cb[11];cb[4];cb[1];cb[8]])
plot!(FR,MAG,label=:false, lw=2, title="Surge Acceleration")
plot!([1/T₀ₛ;1/T₀ₛ],[1e-6;MAGₚ],lab=L"f^{surge}_0",line=:dash, lw=2)
plot!([1/T₀ₚ;1/T₀ₚ],[1e-6;MAGₚ],lab=L"f^{pitch}_0",line=:dash, lw=2)
plot!([1/fₚ;1/fₚ],[1e-6;MAGₚ],lab=L"f^{JONSWAP}_p",line=:dashdot, lw=2)
plot!(legend=:topright, legendcolumns=1)
plot!(xscale=:log10,xlim=(1e-2,10), minorgrid=true)

figid = joinpath(accfigs,"surge_acc_sp")
savefig(plt_SPX,figid*".svg")
savefig(plt_SPX,figid*".png")

NormSurge = MAG./MAGₚ

# Spectral analysis: Acc Y (Sway)
FR, MAG, phi, df, T, Nfft, H = one_side_asp(parsed_cont[:,2], Time)
MAGₚ = findmax(MAG)[1]

plt_SPY = plot(xlab=L"f~[Hz]", ylab=L"Magnitude~[m/s^2]", palette=[cb[11];cb[4];cb[8]])
plot!(FR,MAG,lab="Sway acceleration",lw=2)
plot!([1/T₀ₛ;1/T₀ₛ],[1e-6;MAGₚ],lab=L"f^{surge}_0",line=:dash, lw=2)
plot!([1/fₚ;1/fₚ],[1e-6;MAGₚ],lab=L"f_p~(JONSWAP)",line=:dashdot,lw=2)
plot!(xscale=:log10)
plot!(xlim=(1e-2,10), minorgrid=true)

# Spectral analysis: Acc Z (Heave)
FR, MAG, phi, df, T, Nfft, H = one_side_asp(parsed_cont[:,3].-mean(parsed_cont[:,3]), Time)
MAGₚ = findmax(MAG)[1]

plt_SPZ = plot(xlab=L"f~[Hz]", ylab=L"Magnitude~[m/s^2]", palette=[cb[11];cb[4];cb[8]])
plot!(FR,MAG,label=:false,title="Heave Acceleration", lw=2)
plot!([1/T₀ₕ;1/T₀ₕ],[1e-6;MAGₚ],lab=L"f^{heave}_0",line=:dash, lw=2)
plot!([fₚ;fₚ],[1e-6;MAGₚ],lab=L"f^{JONSWAP}_p",line=:dashdot, lw=2)
plot!(legend=:topright, legendcolumns=1)
plot!(xscale=:log10,xlim=(1e-2,10), minorgrid=true)

figid = joinpath(accfigs,"heave_acc_sp")
savefig(plt_SPZ,figid*".svg")
savefig(plt_SPZ,figid*".png")

NormHeave = MAG./MAGₚ

# Spectral analysis: Ang. Vel Y (Pitch)
FR, MAG, phi, df, T, Nfft, H = one_side_asp(parsed_cont[:,5], Time)
MAGₚ = findmax(MAG)[1]

plt_SPP = plot(xlab=L"f~[Hz]", ylab=L"Magnitude~[rad/s]", palette=[cb[11];cb[4];cb[1];cb[8]])
plot!(FR,MAG,label=:false, title="Pitch Angular Velocity", lw=2)
plot!([1/T₀ₛ;1/T₀ₛ],[1e-6;MAGₚ],lab=L"f^{surge}_0",line=:dash, lw=2)
plot!([1/T₀ₚ;1/T₀ₚ],[1e-6;MAGₚ],lab=L"f^{pitch}_0",line=:dash, lw=2)
plot!([1/fₚ;1/fₚ],[1e-6;MAGₚ],lab=L"f^{JONSWAP}_p",line=:dashdot, lw=2)
plot!(legend=:topright, legendcolumns=1)
plot!(xscale=:log10,xlim=(1e-2,10), minorgrid=true)

figid = joinpath(accfigs,"pitch_rate_sp")
savefig(plt_SPP,figid*".svg")
savefig(plt_SPP,figid*".png")

NormPitch = MAG./MAGₚ

plt_all = plot(xlab="f [Hz]", ylab="Norm. Magnitude [-]", palette=[cb[11];cb[4];cb[8]])
plot!(FR,NormSurge,label="Surge Acceleration")
plot!(FR,NormHeave,label="Heave Acceleration")
plot!(FR,NormPitch,label="Pitch Angular Velocity",line=:dash)
plot!(legend=:topright, legendcolumns=1)
plot!(xscale=:log10,xlim=(1e-2,10), minorgrid=true)

# plt_abc = plot(plt_SPX, plt_SPZ, plt_SPP, layout = @layout [a ; b ; c])

display(plt_SPX)
display(plt_SPY)
display(plt_SPZ)
display(plt_SPP)
display(plt_all)
# display(plt_abc)

###########################################################################
## Cameras
if fcam
    # 1
    plt_c1 = plot(xlab = "t [s]", ylab="[m]", title = "Surge", palette=cb[1:NoMark])
    for i ∈ 1:NoMark
        plot!(t, x[:,i], lw=1, line=:dot, lab = "Marker $i")
    end
    plot!(t, x₀, lw=2, lab="CG")
    plot!(legendcolumns=2, legend=:topright)
    display(plt_c1)

    figid = joinpath(camfigs,"surge_disp")
    savefig(plt_c1,figid*".svg")
    savefig(plt_c1,figid*".png")

    # 2
    plt_c2 = plot(xlab = "t [s]",ylab="[m]", title = "Sway", palette=cb[1:NoMark])
    for i ∈ 1:NoMark
        plot!(t, y[:,i], lw=1, line=:dot, lab = "Marker $i")
    end
    plot!(t, y₀, lw=2, lab="CG")
    plot!(legendcolumns=2, legend=:topright)
    display(plt_c2)

    figid = joinpath(camfigs,"sway_disp")
    savefig(plt_c2,figid*".svg")
    savefig(plt_c2,figid*".png")

    # 3
    plt_c3 = plot(xlab = "t [s]",ylab="[m]", title = "Heave", palette=cb[1:NoMark])
    for i ∈ 1:NoMark
        plot!(t, z[:,i], lw=1, line=:dot, lab = "Marker $i")
    end
    plot!(t, z₀, lw=2, lab="CG")
    plot!(legendcolumns=2, legend=:topright)
    display(plt_c3)

    figid = joinpath(camfigs,"heave_disp")
    savefig(plt_c3,figid*".svg")
    savefig(plt_c3,figid*".png")

    # 4
    plt_cp = plot(xlab = "t [s]",ylab="[°]", title = "Pitch", palette=cb[1:NoMark])
    for i ∈ 1:NoMark
        plot!(t, B.θ[:,i]*180/π, lw=1, line=:dot, lab = "Marker $i")
    end
    plot!(t, θₐ*180/π, lw=2, lab="Accelerometer")
    plot!(legendcolumns=2, legend=:topright)
    display(plt_cp)

    figid = joinpath(camfigs,"pitch_ang")
    savefig(plt_cp,figid*".svg")
    savefig(plt_cp,figid*".png")

    # 5
    plt_c6 = plot(xlab = "x [m]", ylab = "z [m]", title = "X-Z plane", palette=cb[1:NoMark])
    for i ∈ 1:NoMark
        plot!(X[:,i], Z[:,i], lw=1, line=:dot, lab = "Marker $i")
    end
    plot!(CG.x, CG.z, lw=2, lab="CG")
    plot!(legendcolumns=2, legend=:topright)
    display(plt_c6)

    figid = joinpath(camfigs,"X-Z_plane_pos")
    savefig(plt_c6,figid*".svg")
    savefig(plt_c6,figid*".png")

    # 6
    plt_c4 = plot(xlab = "x [m]", ylab = "y [m]", title = "X-Y plane (top view)", palette=cb[1:NoMark])
    for i ∈ 1:NoMark
        plot!(X[:,i], Y[:,i], lw=1, line=:dot, lab = "Marker $i")
    end
    plot!(CG.x, CG.y, lw=2, lab="CG")
    plot!(legendcolumns=2, legend=:topright)
    display(plt_c4)

    figid = joinpath(camfigs,"X-Y_plane_pos")
    savefig(plt_c4,figid*".svg")
    savefig(plt_c4,figid*".png")

    # 7
    plt_c5 = plot(xlab = "y [m]", ylab = "z [m]", title = "Y-Z plane (side view)", palette=cb[1:NoMark])
    for i ∈ 1:NoMark
        plot!(Y[:,i], Z[:,i], lw=1, line=:dot, lab = "Marker $i")
    end
    plot!(CG.y, CG.z, lw=2, lab="CG")
    plot!(legendcolumns=2, legend=:topright)
    display(plt_c5)

    figid = joinpath(camfigs,"Y-Z_plane_pos")
    savefig(plt_c5,figid*".svg")
    savefig(plt_c5,figid*".png")

    # 8
    plt_xzinit = plot(xlab = "x [m]", ylab = "z [m]", title = "X-Z plane (front view)", aspect_ratio=:equal)
    plot!(X[1,:], Z[1,:], seriestype=:scatter, ms=10, mc=cb[1:NoMark], legend=:false)
    for i ∈ 1:NoMark
        annotate!(X[1,i], Z[1,i], text("$i", :center, 10))
    end
    plot!([xG0], [zG0], seriestype=:scatter, ms=10, mc=cb[4], legend=:false)
    annotate!([xG0], [zG0], text("G", :center, 10))
    display(plt_xzinit)

    figid = joinpath(camfigs,"X-Z_init_pos")
    savefig(plt_xzinit,figid*".svg")
    savefig(plt_xzinit,figid*".png")

    # 9
    plt_xyinit = plot(xlab = "x [m]", ylab = "y [m]", title = "X-Y plane (top view)", aspect_ratio=:equal)
    plot!(xᵣ,yᵣ,lw=2,legend=:false)
    plot!(X[1,:], Y[1,:], seriestype=:scatter, ms=10, mc=cb[1:NoMark], legend=:false)
    for i ∈ 1:NoMark
        annotate!(X[1,i], Y[1,i], text("$i", :center, 10))
    end
    plot!([xG0], [yG0], seriestype=:scatter, ms=10, mc=cb[4], legend=:false)
    annotate!([xG0], [yG0], text("G", :center, 10))
    display(plt_xyinit)

    figid = joinpath(camfigs,"Y-X_init_pos")
    savefig(plt_xyinit,figid*".svg")
    savefig(plt_xyinit,figid*".png")

    # 10
    plt_yzinit = plot(xlab = "y [m]", ylab = "z [m]", title = "Y-Z plane (side view)", aspect_ratio=:equal)
    plot!(Y[1,:], Z[1,:], seriestype=:scatter, ms=10, mc=cb[1:NoMark], legend=:false)
    for i ∈ 1:NoMark
        annotate!(Y[1,i], Z[1,i], text("$i", :center, 10))
    end
    plot!([yG0], [zG0], seriestype=:scatter, ms=10, mc=cb[4], legend=:false)
    annotate!([yG0], [zG0], text("G", :center, 10))
    display(plt_yzinit)

    figid = joinpath(camfigs,"Y-Z_init_pos")
    savefig(plt_yzinit,figid*".svg")
    savefig(plt_yzinit,figid*".png")

    # 11
    plt_3Dinit = plot(xlab = "x [m]", ylab = "y [m]", zlab = "z [mm]", title = "Setup", aspect_ratio=:equal)
    plot!(xᵣ,yᵣ, zᵣ, lw=2,legend=:false)
    plot!(X[1,:], Y[1,:], Z[1,:], seriestype=:scatter, ms=10, mc=cb[1:NoMark], legend=:false)
    for i ∈ 1:NoMark
        annotate!(X[1,i], Y[1,i], Z[1,i], text("$i", :center, 10))
    end
    plot!([xG0], [yG0], [zG0], seriestype=:scatter, ms=10, mc=cb[4], legend=:false)
    annotate!([xG0], [yG0], [zG0], text("G", :center, 10))
    display(plt_3Dinit)

    figid = joinpath(camfigs,"3D_init_pos")
    savefig(plt_3Dinit,figid*".svg")
    savefig(plt_3Dinit,figid*".png")

    # 12
    plt_SPc1 = plot(xlab=L"f~[Hz]", ylab=L"Magnitude~[m]", palette=[cb[11];cb[4];cb[1];cb[6];cb[8]])
    plot!(FRx,MAGx,label=:false, title="Surge", lw=2)
    plot!([1/T₀ₛ;1/T₀ₛ],[1e-6;MXₚ],lab=L"f^{surge}_0",line=:dash, lw=2)
    plot!([1/T₀ₚ;1/T₀ₚ],[1e-6;MXₚ],lab=L"f^{pitch}_0",line=:dash, lw=2)
    plot!([1/T₀ₕ;1/T₀ₕ],[1e-6;MXₚ],lab=L"f^{heave}_0",line=:dash, lw=2)
    plot!([1/fₚ;1/fₚ],[1e-6;MXₚ],lab=L"f^{JONSWAP}_p",line=:dashdot, lw=2)
    plot!(xscale=:log10,xlim=(1e-2,10), minorgrid=true)
    display(plt_SPc1)

    figid = joinpath(camfigs,"surge_frf")
    savefig(plt_SPc1,figid*".svg")
    savefig(plt_SPc1,figid*".png")

    # 13
    plt_SPc2 = plot(xlab=L"f~[Hz]", ylab=L"Magnitude~[m]", palette=[cb[11];cb[4];cb[1];cb[6];cb[8]])
    plot!(FRz,MAGz,label=:false, title="Heave", lw=2)
    plot!([1/T₀ₛ;1/T₀ₛ],[1e-6;MZₚ],lab=L"f^{surge}_0",line=:dash, lw=2)
    plot!([1/T₀ₚ;1/T₀ₚ],[1e-6;MZₚ],lab=L"f^{pitch}_0",line=:dash, lw=2)
    plot!([1/T₀ₕ;1/T₀ₕ],[1e-6;MZₚ],lab=L"f^{heave}_0",line=:dash, lw=2)
    plot!([1/fₚ;1/fₚ],[1e-6;MZₚ],lab=L"f^{JONSWAP}_p",line=:dashdot, lw=2)
    plot!(xscale=:log10,xlim=(1e-2,10), minorgrid=true)
    display(plt_SPc2)

    figid = joinpath(camfigs,"heave_frf")
    savefig(plt_SPc2,figid*".svg")
    savefig(plt_SPc2,figid*".png")

    # 14
    plt_SPc3 = plot(xlab=L"f~[Hz]", ylab=L"Angle~[°]", palette=[cb[11];cb[4];cb[1];cb[6];cb[8]])
    plot!(FRθ,MAGθ,label=:false, title="Pitch", lw=2)
    plot!([1/T₀ₛ;1/T₀ₛ],[1e-6;MΘₚ],lab=L"f^{surge}_0",line=:dash, lw=2)
    plot!([1/T₀ₚ;1/T₀ₚ],[1e-6;MΘₚ],lab=L"f^{pitch}_0",line=:dash, lw=2)
    plot!([1/T₀ₕ;1/T₀ₕ],[1e-6;MΘₚ],lab=L"f^{heave}_0",line=:dash, lw=2)
    plot!([1/fₚ;1/fₚ],[1e-6;MΘₚ],lab=L"f^{JONSWAP}_p",line=:dashdot, lw=2)
    plot!(xscale=:log10,xlim=(1e-2,10), minorgrid=true)
    display(plt_SPc3)

    figid = joinpath(camfigs,"pitch_frf")
    savefig(plt_SPc3,figid*".svg")
    savefig(plt_SPc3,figid*".png")

    # 15
    pltRAOs = plot(xlab=L"f~[Hz]", ylab=L"RAO~[-]", palette=[cb[8];cb[11];cb[4]])
    plot!(FRx,MAGx./MAGx[2], lab="Surge", lw=2)
    plot!(FRz,MAGz./MAGz[2], lab="Heave", lw=2)
    plot!(FRθ,MAGθ./MAGθ[2], lab="Pitch", lw=2)
    plot!([fˣ;fˣ],[1e-6;MXₚ./MAGx[2]],lab=L"f^{surge}_0",line=:dash, lw=1)
    plot!([fᶻ;fᶻ],[1e-6;MXₚ./MAGx[2]],lab=L"f^{pitch}_0",line=:dash, lw=1)
    plot!([fθ;fθ],[1e-6;MXₚ./MAGx[2]],lab=L"f^{heave}_0",line=:dash, lw=1)
    plot!(xscale=:log10,xlim=(1e-2,10), minorgrid=true)
    plot!(ylim=(0,MXₚ./MAGx[2]))
    display(pltRAOs)

    figid = joinpath(camfigs,"RAOs_comp")
    savefig(pltRAOs,figid*".svg")
    savefig(pltRAOs,figid*".png")

    # 16
    plt_vel = plot(xlab = "t [s]", ylab = "[m/s]", title="CG velocity", palette=[cb[11];cb[4];cb[8]])
    plot!(t, u₀, lw=2, lab=L"u_0")
    plot!(t, v₀, lw=2, line=:dash, lab=L"v_0")
    plot!(t, w₀, lw=2, line=:dot, lab=L"w_0")
    display(plt_vel)

    plt_acc = plot(xlab = "t [s]", ylab = "[m/s^2]", title="CG acceleration", palette=[cb[11];cb[4];cb[8]])
    plot!(t, u̇₀, lw=2, lab=L"a_x")
    plot!(t, v̇₀, lw=2, line=:dash, lab=L"a_y")
    plot!(t, ẇ₀, lw=2, line=:dot, lab=L"a_z")
    display(plt_acc)

    plt_ω = plot(xlab="[s]", ylab="[rad/s]", title="Angular velocities", palette=[cb[11];cb[4];cb[8]])
    plot!(t,ωˣ, lw=2, line=:dash, lab=L"\omega_x")
    plot!(t,ωʸ, lw=2, lab=L"\omega_y")
    plot!(t,ωᶻ, lw=2, line=:dot, lab=L"\omega_z")
    display(plt_ω)

    if fprb
        mid = 1;
        plt_compX = plot(ylab="[m]", palette=[cb[11];cb[1]])
        plot!(t, x[:,mid], lw=2, lab = "Surge (M #$mid)")
        plot!(tₚᵣ,Zₚᵣ,lw=2,lab="Probe #4")
        plot!(xlim=(0,tₚᵣ[end]))

        plt_compZΘ = plot(palette=[cb[1];cb[8];cb[4]])
        plot!(tₚᵣ,Zₚᵣ,lw=2,lab="Probe #4")
        plot!(t, z[:,mid], lw=2, lab = "Heave [m] (M #$mid)")
        plot!(t, B.θ[:,mid]*180/π, lw=2, lab = "Pitch [°] (M #$mid)")
        plot!(xlim=(0,tₚᵣ[end]))

        ## Probe 4
        plt_prb4 = plot(xlab="t [s]", ylab="η [m]", title = "Surface elevation", palette=[cb[1]])
        plot!(tₚᵣ,Zₚᵣ,lw=2,lab="Probe #4")
        display(plt_prb4)

        figid = joinpath(prbfigs,"clb_pr4_figs",dname)
        savefig(plt_prb4,figid*".svg")
        savefig(plt_prb4,figid*".png")

        plt_resp = plot(plt_compX, plt_compZΘ, layout = @layout [a; b])
        display(plt_resp)

        figid = joinpath(camfigs,"elev_resp_comp")
        savefig(plt_resp,figid*".svg")
        savefig(plt_resp,figid*".png")
    end
    display(plt_fairten)
end
function parse_fxw(fid, hdl)
    # Parse fixed width data from a .txt file in the library folder
    # fid: file path,   
    # skiphead = 0,1: The file has a header that must be skipped (No,Yes)


    open(fid, "r")
    if hdl == 0
        content = readdlm(fid, '\n')
    else
        content = readdlm(fid, skipstart=hdl, '\n')
        println("$hdl header lines have been skipped.")
    end

    NoRows = length(content)
    NoCols = length(split(content[1]))

    Parsed_content = zeros(Float64,NoRows, NoCols)

    for i ∈ 1:NoRows
        for j ∈ 1:NoCols
            tmp = split(content[i])
            Parsed_content[i,j] = parse(Float64, tmp[j])
        end
    end

    return Parsed_content
end

function parse_fxw_pf(fid, w_new, hdl)
    ## Parse fixed width data from a text file
    # fid: file path,   
    # w_new=0,1: write parsed data to new file (No,Yes),
    # hdl: Header lines to be skipped

    fid_new  = fid * "_new" # Path to new file

    open(fid, "r")
    if hdl == 0
        content = readdlm(fid, '\n')
    else
        content = readdlm(fid, skipstart=hdl, '\n')
        println("$hdl header lines have been skipped.")
    end

    NoRows = length(content)
    NoCols = length(split(content[1]))

    Parsed_content = zeros(Float64,NoRows, NoCols)

    for i ∈ 1:NoRows
        for j ∈ 1:NoCols
            tmp = split(content[i])
            Parsed_content[i,j] = parse(Float64, tmp[j])
        end
    end

    if w_new == 1
        writedlm(fid_new, Parsed_content, '\t')
    end

    return Parsed_content
end

function parse_fixed_width(fid, w_new)
    # Parse fixed width data from a .txt file in the library folder
    # fid: file path,   
    # w_new=0,1: write parsed data to new file (No,Yes),
    # skiphead = 0,1: The file has a header that must be skipped (No,Yes)

    fid_new  = fid[1:end-4] * "_new.txt" # Path to new file

    open("library/"*fid, "r")
    content = readdlm("library/"*fid, '\n')

    NoRows = length(content)
    NoCols = length(split(content[1]))

    Parsed_content = zeros(Float64,NoRows, NoCols)

    for i ∈ 1:NoRows
        for j ∈ 1:NoCols
            tmp = split(content[i])
            Parsed_content[i,j] = parse(Float64, tmp[j])
        end
    end

    if w_new == 1
        writedlm("library/"*fid_new, Parsed_content, '\t')
    end

    return Parsed_content
end

function read_ow3d_inp(fid)
    open("library/"*fid, "r")
    cont = readdlm("library/"*fid, '\n')

    # 1: Header
    Head = cont[1]

    len = length(cont)
    AllValues = zeros(len-1,13).*NaN
    AllStrings = fill(" ", len-1, 13)

    for i ∈ 2:len
        if length(cont[i]) == 1
            AllValues[i-1,1] = cont[i]
        else
            tmp_str_vec = split(cont[i]) # Temporary string vector
            Ltmp = length(tmp_str_vec) # Temp length
        
            for j ∈ 1:Ltmp
                if isnothing(tryparse(Float64,tmp_str_vec[j]))
                    AllStrings[i-1,j] = tmp_str_vec[j]
                else
                    AllValues[i-1,j] = parse(Float64,tmp_str_vec[j])
                    AllStrings[i-1,j] = "-"
                end
            end
        end
    end

    curr_time = Dates.format(now(), "HH:MM:SS")
    print("--------------------------------------------------------\n")
    print("OW3D Input Data: HH:MM:SS = ", curr_time, "\n")

    # 2: Initial Conditions
    InitConds = AllValues[1,:]
    IC, IncWT = [Int(InitConds[i]) for i ∈ 1:2]
    acc_tol = InitConds[3]

    if IC == -1
        print("-> Initial Conditions: Read from 'OceanWave3D.init' , Not implemented curvilinear!\n")
    elseif IC == 0
        print("-> Initial Conditions: Determined by PressureTerm\n")
    else
        print("-> Initial Conditions: Specific case, non-standard definition\n")
    end

    if IncWT == 3
        print("-> Wave Type = Wavemaker \n")
        Tₚ, Hₛ, d, γ, max_kd = 8, 3, 100, 3.3, 157
    end

    # 3: Computational domain and resolution parameters
    Domain = AllValues[2,1:6]
    Lx, Ly, Lz = [Domain[i] for i ∈ 1:3]
    Nx, Ny, Nz = [Int(Domain[i]) for i ∈ 4:6]

    if Ny > 1
        print("-> Geometry = 3D \n")
    else
        print("-> Geometry = 2D \n")
        print("-> Computational Domain: ", "Lₓ=", Lx, " [m], ", "Nₓ=", Nx,", ", "Nz=", Nz, "\n")
    end

    # 4: Finite difference - Preconditioner

    # 5: Time parameters
    Time = AllValues[4,:]
    Nt, dt = Int(Time[1]), Time[2]
    print("-> Time parameters: ", "Nₜ=", Nt, ", ", "Δt=", dt, " [s]\n")
    
    # 6: Global constants
    g = AllValues[5,1]; ρw = AllValues[5,2]

    # 7: Solver parameters
    # 8: Stream function solution parameters
    # 9: Data storage info
    # 10: Domain & Time resolution start/stop/stride
    # 11: Relaxation zone setup - Output file tmin,tmax

    # 12: Mode, applied free-surface pressure
    NonLin, PressTerm = AllValues[11,1], AllValues[11,2]
    if NonLin == 1
        print("-> Mode: Fully nonlinear model is employed.\n")
    elseif NonLin == 0
        print("-> Mode: Linear model is employed.\n")
    else
        print("ERROR: Could not identify model Mode! Value missing or wrongly defined. \n")
    end

    if PressTerm == 0
        print("-> Free surface pressure: No term is being applied\n")
    elseif PressTerm == 1
        print("-> Free surface pressure: 2D Gaussian\n")
    elseif PressTerm == 2
        print("-> Free surface pressure: 3D Gaussian\n")
    elseif PressTerm == 3
        print("-> Free surface pressure: 3D tanh\n")
    elseif PressTerm == 4
        print("-> Free surface pressure: 3D cos²\n")
    else
        print("WARNING: No field defined. \n")
    end

    # 13: SG-FILTERING
    # 14: Relaxation zones

    # 15: Wave Generation Zone
    WaveGenZone = AllValues[14,:]
    xbeg_WGZ, xend_WGZ, ybeg_WGZ, yend_WGZ = [WaveGenZone[i] for i ∈ 1:4]
    print("-> Wave Generation Zone: ", "x = ", xbeg_WGZ, "-", xend_WGZ, " [m]\n")

    # 16: Wave Absorption Zone
    WaveAbsZone = AllValues[15,:]
    xbeg_WAZ, xend_WAZ, ybeg_WAZ, yend_WAZ = [WaveAbsZone[i] for i ∈ 1:4]
    print("-> Wave Absorption Zone: ", "x = ", xbeg_WAZ, "-", xend_WAZ, " [m]\n")

    # 17: Damping Pressure Zones
    PressDamping, NoPDZns = AllValues[16,1], AllValues[16,2]

    if PressDamping == 0
        print("-> Pressure damping: OFF\n")
    elseif PressDamping == 1
        print("-> Pressure damping: ON, No of PD zones=", NoPDZns, "\n")
    end

    # 18: Damping pressure zone parameters

    # 19: Curvelinear coordinates
    Curvelin = AllValues[18,1]

    if Curvelin == 0
        print("-> Standard Cartesian model employed. \n")
    elseif Curvelin == 1
        print("-> Curvelinear model employed. \n")
    else
        print("ERROR: Could not identify if Curvelinear is ON or OFF! Value missing or wrongly defined. \n")
    end

    # 20: Wave type parameters
    WaveType = AllValues[end,:]
    WT_id =WaveType[1] 

    if WT_id == 1
        print("-> Wave Type = Linear Wave Theory \n")
    elseif WT_id == 2
        print("-> Wave Type = Stream Function \n")
    elseif WT_id == 3
        Tₚ, Hₛ, d, max_kd, seed1, seed2, xGW, yGW, γ = [WaveType[i] for i ∈ 2:10]
        print("-> Wave Type = JONSWAP \n",
            "Tₚ=",Tₚ," [s], ","Hₛ=",Hₛ," [m], ", "d=", d, " [m], ", "max(kd)=", max_kd, ", ", "γ=", γ, "\n",
            "The generated wave is centered at (x,y)=", "(", xGW, ",", yGW, ")", "[m]\n")
    elseif WT_id == 4
        Tₚ, Hₛ, d, max_kd, seed1, seed2, xGW, yGW = [WaveType[i] for i ∈ 2:9]
        γ = 1
        print("-> Wave Type = Pierson-Moskowitz \n",
        "Tₚ=",Tₚ," [s]\t","Hₛ=",Hₛ," [m]\t", "d=", d, " [m]\t", "max(kd)=", max_kd, "\n",
        "The generated wave is centered at (x,y)=", "(", xGW, ",", yGW, ")", " [m]\n")
    elseif WT_id == 5
        print("-> Wave Type = Wave File \n")
    elseif WT_id == 6
        print("-> Wave Type = Wavemaker \n")
        Tₚ, Hₛ, d, γ = 8, 3, 100, 3.3
    else
        print("ERROR: Could not identify Wave Type! Value missing or wrongly defined. \n")
    end

    print("--------------------------------------------------------\n")

    return Lx, Nx, Nz, Tₚ, Hₛ, γ, d, max_kd, g, ρw
end
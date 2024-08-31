function paths(pdir, run_id, phi_id, prb_id, flong)
    # Definition of necessary paths and formulation of requested directories and file paths
    # Inputs:
    #   pdir::String = Parent directory in library folder ("SE", "WM", "rand_seas" etc.)
    #   run_id::Int = Run identifier in each case 
    #   phi_id::Int = Phase shift for decomposition: 0, 1, 2 , 3 == rand(), -π/2, -π, -3π/2 
    #   prb_id::Int = ID No of probe from OW3D simulations (see var. probes below)
    #   flong::Boolean = Flag for long OW3D simulations

    # Outputs:
    #   phipath::String = Path to phase shifted realization for specified run under selected case
    #   Dpath::String = Path to Decomposition folder
    #   DecFigs::String = Path to Decomposition figures
    #   DecEvs::String = Path to events processed from Decomposition
    #   postOFpath::String = Path OpenFast postprocessing folder 
    #   postOW3Dpath::String = Path to OW3D postprocessing folder
    #############################################################################################
    ## Definition of absolute parent directory paths containg necessary data
    HOME::String = "/home/andreasp"                 # Path to home directory
    HardDrive::String = "/run/media/data_storage"   # Path to secondary hard drive
    ProjPath::String = pwd()                        # Project path

    libpath::String = ProjPath*"/library" # Path to project's library folder
    OW3Dpath::String = HOME*"/OW3D_Out/WAVEGEM"   # Path to OW3D simulation output folder
    OFASTpath::String = HardDrive*"/ONISILOS/OFAST"   # Path to OW3D simulation output folder

    probes = ("x0", "x100", "x700", "x1000") # OW3D probes

    #############################################################################################
    ## Automate library directory
    case_id = js_case(Hₛ, Tₚ, γ, flong) # Case of JONSWAP realisation (set of Tp, Hs, γ)
    head_case = "# JS_$(Hₛ)m_$(Tₚ)s_ϕ$phi_id"
    OW3Dprb = probes[prb_id]

    runstr::String = "/$(pdir)/$(case_id)/$(run_id)"
    casedir::String = libpath * "/$(pdir)/$(case_id)/"
    rundir::String = libpath * runstr
    phipath::String = libpath * runstr * "/$(phi_id)/"
    OW3Dcdir::String = OW3Dpath * "/$(pdir)/$(case_id)/"
    OW3Drdir::String = OW3Dpath * runstr
    OW3Dphipath::String = OW3Dpath * runstr * "/$(phi_id)/"
    Decpath::String = libpath * "/$(pdir)/$(case_id)/$(run_id)/Decomposition/"
    DecFigs::String = Decpath * "Figures/"
    DecEvs::String = Decpath * "events/"
    postOFpath::String = phipath * "postOFAST/"
    postOW3Dpath::String = phipath * "postOW3D/$(OW3Dprb)/"

    # println("Paths:")
    # println("Processing the case: $case_id")
    # println("Case specifics: $head_case")
    # println("Run folder: $runstr")
    # println("Path to OW3D input folder: $OW3Dinp")
    # println("Case path in library folder: $phipath")
    # println("Case path in OW3D folder: $postOW3Dpath")
    # println("Case path in OFAST folder: $postOFpath")

    return runstr, case_id, casedir, rundir, phipath, OW3Dcdir, OW3Drdir, OW3Dphipath, 
            Decpath, DecEvs, DecFigs, postOFpath, postOW3Dpath, OFASTpath
end

function make_dirs(levels, folders...)
    # Creates necessary new directories
    if levels == 3
        parent = folders[1]
        child = folders[2]
        gchild = folders[3]

        if !isdir(parent)
            mkdir(parent)
            mkdir(child)
            mkdir(gchild)
        elseif !isdir(child)
            mkdir(child)
            mkdir(gchild)
        elseif !isdir(gchild)
            mkdir(gchild)
        else
            println("WARNING: Folder already exists. Contained files may be overwritten.")
            println("WARNING: Are you sure you want to proceed? (y/n)")
            uinp = readline()
            if uinp == "n"
                error("Script execution aborted!")
            end
        end
    elseif levels == 2
        parent = folders[1]
        child = folders[2]

        if !isdir(parent)
            mkdir(parent)
            mkdir(child)
        elseif !isdir(child)
            mkdir(child)
        else
            println("WARNING: Folder already exists. Contained files may be overwritten.")
            println("WARNING: Are you sure you want to proceed? (y/n)")
            uinp = readline()
            if uinp == "n"
                error("Script execution aborted!")
            end
        end
    elseif levels == 1
        parent = folders[1]

        if !isdir(parent)
            mkdir(parent)
        else
            println("WARNING: Folder already exists. Contained files may be overwritten.")
            println("WARNING: Are you sure you want to proceed? (y/n)")
            uinp = readline()
            if uinp == "n"
                error("Script execution aborted!")
            end
        end
    else
        println("Not applicable level of folders has been provided. Nothing will happen.")
    end
    
    return nothing
end

function js_case(Hₛ, Tₚ, γ, flag)
    Tpγ_id = ""
    Hs_id = ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')

    Hs_sch = findall( x -> x == Hₛ, collect(range(1,8)))

    if (length(Hs_sch) == 0) 
        println("Custom case of Hₛ.")
        case_id = "JS_$(Hₛ)m_$(Tₚ)s"
    else
        if Tₚ == 8 && γ == 3.3
            Tpγ_id = "A"
        elseif Tₚ == 9 && γ == 3.3
            Tpγ_id = "B"
        elseif Tₚ == 10 && γ == 3.3
            Tpγ_id = "C"
        elseif Tₚ == 8 && γ == 2.0
            Tpγ_id = "D"                    
        elseif Tₚ == 8 && γ == 5.0
            Tpγ_id = "E"
        else
            println("Custom case of JONSWAP realisation.")
        end
        case_id = Tpγ_id*Hs_id[Int(Hₛ)]

        if flag
            case_id = "L"*case_id
        end
    end

    return case_id
end
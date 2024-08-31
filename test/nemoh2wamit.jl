using Printf
using DelimitedFiles

include("text_process.jl")

# fid1 = "/home/andreasp/NEMOH/Run_00/mesh/KH.dat"
fid1 = "/home/andreasp/NEMOH/Mike/Mechanics/Kh.dat"
fid2 = "/home/andreasp/NEMOH/Run_00/results/CM.dat"
fid3 = "/home/andreasp/NEMOH/Run_00/results/CA.dat"

const ρ, g::Float64 = 1025.0, 9.81  # Sea water density [kg/m³], Accel. of gravity [m/s²]

KH = kh2hst(fid1)
# hcoef2wam1(fid2)
# hcoef2wam1(fid3)

function hcoef2wam1(fid)
    open(fid, "r")
    cont = readdlm(fid, '\n')

    # 1: Header
    Head = cont[1]
    Nom = parse(Int, split(Head)[end])

end

function kh2hst(fid)
    w_new = 0
    hdl = 0

    cont = parse_fxw_pf(fid, w_new, hdl)
    cont = cont./(ρ*g)

    NoRows = size(cont)[1]
    NoCols = size(cont)[2]

    KH = (x->(@sprintf "%.7e" x)).(cont)

    if NoRows != NoCols || NoRows != 6 || NoCols != 6
        error("The hydrostatic matrix must be 6x6.")
    else
        fid_new  = fid[1:end-3] * "hst" # Path to new file
        open(fid_new, "w")
        OutTpl = Vector{Tuple}(undef, NoRows^2)
        cnt = 1

        for i ∈ 1:NoRows
            for j ∈ 1:NoCols
                OutTpl[cnt] = (i, j, KH[i,j]) 
                cnt += 1
            end
        end
    end

    writedlm(fid_new, OutTpl, '\t')

    return KH
end
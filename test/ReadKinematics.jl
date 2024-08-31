# This script reads the unformatted binary kinematics output file from the OceanWave3D code.
using Plots, LaTeXStrings, DelimitedFiles, Dates, LinearAlgebra
using FFTW, DSP, Statistics
using Printf
using SparseArrays

include("directories.jl")

# INPUT
Nbits::Int8 = 32  # Set Nbits=32 or 64
## Significant wave height [m], Peak period [s], peakedness [-], Cut-off frequency [Hz]
Hₛ, Tₚ, γ, fcut::Float64 = 6.0, 8.0, 3.3, 0.625  # JONSWAP parameters 
pdir::String = "SE" # Parent directory in library folder ("SE", "WM", "rand_seas" etc.)
run_id::Int8 = 1    # Run identifier in each case (1,2,3..)
phi_id::Int8 = 0    # Phase shift for decomposition: 0, 1, 2 , 3 == rand(), -π/2, -π, -3π/2
prb_id::Int8 = 4    # ID No of probe from OW3D simulations
flong = Bool(1) # Long simulation? Set 0: false or 1: true

runstr, case_id, casedir, rundir, phipath, OW3Dcdir, OW3Drdir, OW3Dphipath, DecPath,
            _, _, _, postOW3Dpath, _ = paths(pdir, run_id, phi_id, prb_id, flong)

# Paths
# dirp = "SE/"
# casep = "LAE/"
# runp = "2/"
# phsp = "0"

# jullib = "/home/andreasp/WAVEGEM/library/"
# post = "/postOW3D/"
# decomp = "Decomposition/"

# run_folder = joinpath(dirp, casep, runp, phsp)  # Folder of OCW3D run results
# save_fld = joinpath(jullib, run_folder, post)  # Path to Julia lib folder
# save_fld2 = joinpath(jullib, dirp, casep, runp, decomp)  # Path to Decomposition

WriteOut = "yes"  # "yes","no": Save variables to text files
Pressure = "yes"  # "yes","no": Compute time-derivatives and pressures
PlotVar = "yes"  # "yes","no": Plot data

# INPUT CHECKS
if !isdefined(Main, :Nbits) || !isdefined(Main, :prb_id)
    error("Set 'Nbits'=32 or 64 and 'prb_id'= kinematics file number to read")
end

# Open the file
fname = prb_id < 10 ? joinpath(OW3Dphipath, "Kinematics0$(prb_id).bin") : joinpath(OW3Dphipath, "Kinematics$(prb_id).bin")
fid1 = open(fname, "r")

# Check for problems
if fid1 === nothing
    error("Could not open file $fname, returned fid1 = - 1")
end

# Choose number of bits
if Nbits == 32
    int_nbit = Int32
    njunkread = 2
elseif Nbits == 64
    int_nbit = Int64
    njunkread = 2
else
    println("Illegal value for Nbits, Nbits = $Nbits")
end

# SCRIPT
# Read the data from the file
junk = read(fid1, int_nbit)
xbeg = read(fid1, Int32)
xend = read(fid1, Int32)
xstride = read(fid1, Int32)
ybeg = read(fid1, Int32)
yend = read(fid1, Int32)
ystride = read(fid1, Int32)
tbeg = read(fid1, Int32)
tend = read(fid1, Int32)
tstride = read(fid1, Int32)
dt = read(fid1, Float64)
nz = read(fid1, Int32)

junk = read(fid1, int_nbit)
junk = read(fid1, int_nbit)

nx = div((xend - xbeg), xstride) + 1
ny = div((yend - ybeg), ystride) + 1
nt = div((tend - tbeg), tstride) + 1

# A scratch vector for reading the data
tmp = zeros(Float64, nx * ny * max(nz, 5))

# The x-y grid, the depth and bottom gradients for this slice of data
for i ∈ 1:5*nx*ny
    tmp[i] = read(fid1, Float64)
end

junk = read(fid1, int_nbit)
junk = read(fid1, int_nbit)

𝚶⃗ₓ = zeros(Float64, nx, ny)
x, y, h, hx, hy = (𝚶⃗ₓ[:] for _ = 1:5)
x[:] = tmp[1:5:5*nx*ny]
y[:] = tmp[2:5:5*nx*ny]
h[:] = tmp[3:5:5*nx*ny]
hx[:] = tmp[4:5:5*nx*ny]
hy[:] = tmp[5:5:5*nx*ny]

# The sigma coordinate
sigma = zeros(Float64, nz)
for i in 1:nz
    sigma[i] = read(fid1, Float64)
end
junk = read(fid1, int_nbit)
junk = read(fid1, int_nbit)

# Initialize arrays for the solution on this slice
eta = zeros(Float64, nt, nx, ny)
etax = zeros(Float64, nt, nx, ny)
etay = zeros(Float64, nt, nx, ny)
phi = zeros(Float64, nt, nz, nx, ny)
w = zeros(Float64, nt, nz, nx, ny)
u = zeros(Float64, nt, nz, nx, ny)
uz = copy(u)
ux = copy(u)
uy = copy(u)
v = copy(u)
vz = copy(u)
vx = copy(v)
vy = copy(v)
wz = copy(w)
wx = copy(w)
wy = copy(w)
t = [0:nt-1] * dt * tstride  # The time axis

# Read in the solution variables eta, gradeta, phi, u, v, w, dudz, dvdz.
for it in 1:nt
    try
        tmp[1:nx*ny] = read(fid1, Float64, nx*ny)
        eta[it, :, :] = reshape(tmp[1:nx*ny], nx, ny)
        junk = read(fid1, int_nbit, 2)
    catch
        @warn "Read failed at time step $it"
        break
    end

    tmp[1:nx*ny] = read(fid1, Float64, nx*ny)
    etax[it, :, :] = reshape(tmp[1:nx*ny], nx, ny)
    junk = read(fid1, int_nbit, 2)

    tmp[1:nx*ny] = read(fid1, Float64, nx*ny)
    etay[it, :, :] = reshape(tmp[1:nx*ny], nx, ny)
    junk = read(fid1, int_nbit, 2)

    tmp = read(fid1, Float64, nx*ny*nz)
    phi[it, :, :, :] = reshape(tmp, nz, nx, ny)
    junk = read(fid1, int_nbit, 2)

    tmp = read(fid1, Float64, nx*ny*nz)
    u[it, :, :, :] = reshape(tmp, nz, nx, ny)
    junk = read(fid1, int_nbit, 2)

    tmp = read(fid1, Float64, nx*ny*nz)
    v[it, :, :, :] = reshape(tmp, nz, nx, ny)
    junk = read(fid1, int_nbit, 2)

    tmp = read(fid1, Float64, nx*ny*nz)
    w[it, :, :, :] = reshape(tmp, nz, nx, ny)
    junk = read(fid1, int_nbit, 2)

    tmp, count = read(fid1, Float64, nx*ny*nz)
    # Check for an incomplete run.
    if count == 0
        it -= 1
        break
    end
    wz[it, :, :, :] = reshape(tmp, nz, nx, ny)
    junk = read(fid1, int_nbit, 2)

    tmp, count = read(fid1, Float64, nx*ny*nz)
    if count == 0
        it -= 1
        break
    end
    wx[it, :, :, :] = reshape(tmp, nz, nx, ny)
    junk = read(fid1, int_nbit, 2)

    tmp, count = read(fid1, Float64, nx*ny*nz)
    if count == 0
        it -= 1
        break
    end
    wy[it, :, :, :] = reshape(tmp, nz, nx, ny)
    junk = read(fid1, int_nbit, 2)

    tmp, count = read(fid1, Float64, nx*ny*nz)
    if count == 0
        it -= 1
        break
    end
    uz[it, :, :, :] = reshape(tmp, nz, nx, ny)
    junk = read(fid1, int_nbit, 2)

    tmp, count = read(fid1, Float64, nx*ny*nz)
    if count == 0
        it -= 1
        break
    end
    ux[it, :, :, :] = reshape(tmp, nz, nx, ny)
    junk = read(fid1, int_nbit, 2)

    tmp, count = read(fid1, Float64, nx*ny*nz)
    if count == 0
        it -= 1
        break
    end
    uy[it, :, :, :] = reshape(tmp, nz, nx, ny)
    junk = read(fid1, int_nbit, 2)

    tmp, count = read(fid1, Float64, nx*ny*nz)
    if count == 0
        it -= 1
        break
    end
    vz[it, :, :, :] = reshape(tmp, nz, nx, ny)
    junk = read(fid1, int_nbit, 2)

    tmp, count = read(fid1, Float64, nx*ny*nz)
    if count == 0
        it -= 1
        break
    end
    vx[it, :, :, :] = reshape(tmp, nz, nx, ny)
    junk = read(fid1, int_nbit, 2)

    tmp, count = read(fid1, Float64, nx*ny*nz)
    if count == 0
        it -= 1
        break
    end
    vy[it, :, :, :] = reshape(tmp, nz, nx, ny)
    junk = read(fid1, int_nbit, 2)
end

close(fid1)
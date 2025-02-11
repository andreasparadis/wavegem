# Load necessary modules
using Plots, LaTeXStrings, LinearAlgebra, DelimitedFiles
using FFTW, DSP, Statistics, BSplineKit
import ColorSchemes.darkrainbow

gr(fontfamily = "Computer Modern", titlefont = (11, "New Century Schoolbook Bold"))
cb = darkrainbow

# Include necessary scripts for functions
include("text_process.jl"); include("signal_processing.jl")
include("peak_detect.jl");  include("runge_kutta.jl")
include("Gauss_Reg_t.jl");  include("directories.jl")
include("wave_theory.jl");  include("jonswap.jl")

function WaveElevTimeSeriesAtXY_Diff(
    Xcoord::Float64,
    Ycoord::Float64,
    WaveElevSeriesAtXY::Vector{Float64},
    ErrStatLcl::Ref{Int},
    ErrMsgLcl::Ref{String}
)
    # Local variables
    n, m, mu_minus, NStepWave, NStepWave2 = 0, 0, 0, InitInp.NStepWave, InitInp.NStepWave2
    k_n, k_m, L_minus, R_n, R_m, Omega_n, Omega_m, D_minus = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    TmpFreqSeries = Complex{Float64}[0.0 + 0.0im for _ in 1:NStepWave2]
    WaveElevSeriesAtXY .= 0.0

    ErrMsgLcl[] = ""
    ErrStatLcl[] = 0  # Assuming ErrID_None is 0

    for mu_minus in 1:NStepWave2 - 1
        Omega_minus = mu_minus * InitInp.WaveDOmega
        if Omega_minus >= InitInp.WvLowCOffD && Omega_minus <= InitInp.WvHiCOffD
            for m in 1:NStepWave2 - mu_minus
                n = mu_minus + m
                Omega_n = n * InitInp.WaveDOmega
                Omega_m = m * InitInp.WaveDOmega

                k_n = WaveNumber(Omega_n, InitInp.Gravity, InitInp.WtrDpth)
                k_m = WaveNumber(Omega_m, InitInp.Gravity, InitInp.WtrDpth)
                R_n = k_n * tanh(k_n * InitInp.WtrDpth)
                R_m = k_m * tanh(k_m * InitInp.WtrDpth)

                D_minus = TransFuncD_minus(n, m, k_n, k_m, R_n, R_m)

                L_minus = (
                    (D_minus - k_n * k_m *
                     cosd(InitInp.WaveDirArr[n] - InitInp.WaveDirArr[m]) - R_n * R_m) /
                    sqrt(R_n * R_m) + R_n + R_m
                ) / 4.0

                WaveElevxyPrime0 = exp(-1im *
                    ((k_n * cosd(InitInp.WaveDirArr[n]) - k_m * cosd(InitInp.WaveDirArr[m])) * Xcoord +
                     (k_n * sind(InitInp.WaveDirArr[n]) - k_m * sind(InitInp.WaveDirArr[m])) * Ycoord))

                WaveElevC_n = WaveElevC0Norm[n]
                WaveElevC_m = WaveElevC0Norm[m]

                TmpFreqSeries[mu_minus] += 2.0 * WaveElevC_n * conj(WaveElevC_m) * L_minus * WaveElevxyPrime0
            end
        end
    end

    TmpFreqSeries .= TmpFreqSeries ./ 2.0

    # Apply inverse FFT
    ErrStatLcl2 = 0
    ApplyFFT_cx!(WaveElevSeriesAtXY, TmpFreqSeries, FFT_Data, ErrStatLcl2)

    if ErrStatLcl2 != 0
        SetErrStat!(
            ErrStatLcl2,
            "Error occurred while applying the FFT on WaveElevSeriesAtXY.",
            ErrStatLcl,
            ErrMsgLcl,
            "WaveElevTimeSeriesAtXY_Diff"
        )
    end

    WaveElevSeriesAtXY[NStepWave] = WaveElevSeriesAtXY[1]
end

function TransFuncD_minus(n::Int, m::Int, k_n::Float64, k_m::Float64, R_n::Float64, R_m::Float64)::Float64
    # Local variables
    k_nm = 0.0
    SqrtRnMinusRm = 0.0
    Den = 0.0
    Num1 = 0.0
    Num2 = 0.0

    # If n == m, D_minus is set to zero
    if n == m
        return 0.0
    else
        # Compute k_nm
        k_nm = k_nm_minus(n, m, k_n, k_m)

        # Calculate R_nm
        SqrtRnMinusRm = sqrt(R_n) - sqrt(R_m)

        # Calculate the two terms of the numerator
        Num1 = SqrtRnMinusRm * (sqrt(R_m) * (k_n^2 - R_n^2) - sqrt(R_n) * (k_m^2 - R_m^2))
        Num2 = 2 * (SqrtRnMinusRm^2) * (k_n * k_m * cos(deg2rad(InitInp.WaveDirArr[n]) - deg2rad(InitInp.WaveDirArr[m])) + R_n * R_m)

        # Calculate the denominator
        Den = SqrtRnMinusRm^2 - k_nm * tanh(k_nm * InitInp.WtrDpth)

        # Calculate and return TransFuncD_minus
        TransFuncD_minus = (Num1 + Num2) / Den
        return TransFuncD_minus
    end
end

function k_nm_minus(n::Int, m::Int, k_n::Float64, k_m::Float64)::Float64
    # Local variables
    k_nm_minus = 0.0

    if n == m
        # To eliminate any numerical error
        return 0.0
    else
        # Added abs() to handle small negative values
        k_nm_minus = sqrt(abs(k_n^2 + k_m^2 - 2 * k_n * k_m * cos(deg2rad(InitInp.WaveDirArr[n]) - deg2rad(InitInp.WaveDirArr[m]))))
        return k_nm_minus
    end
end

function WaveNumber(Omega::Float64, g::Float64, h::Float64)::Float64
    # Local variables
    C = 0.0
    CC = 0.0
    C2 = 0.0
    A = 0.0
    B = 0.0
    X0 = 0.0
    E2 = 0.0

    # Handle the case when Omega is zero
    if Omega == 0.0
        return 0.0
    else
        # Compute C and CC
        C = (Omega^2 * h) / g
        CC = C^2

        # Compute the initial guess X0
        if C <= 2.0
            X0 = sqrt(C) * (1.0 + C * (0.169 + 0.031 * C))
        else
            E2 = exp(-2.0 * C)
            X0 = C * (1.0 + E2 * (2.0 - 12.0 * E2))
        end

        # Compute the WaveNumber based on the conditions
        if C <= 4.8
            C2 = CC - X0^2
            A = 1.0 / (C - C2)
            B = A * ((0.5 * log((X0 + C) / (X0 - C))) - X0)
            return (X0 - B * C2 * (1.0 + A * B * C * X0)) / h
        else
            return X0 / h
        end
    end
end

function cosh_num_over_cosh_den(k::Float64, h::Float64, z::Float64)::Float64
    """
    Computes the hyperbolic cosine numerator over denominator:
    cosh(k * (z + h)) / cosh(k * h)

    Parameters:
    - k: Wave number (k >= 0) [1/m]
    - h: Water depth (h > 0) [m]
    - z: Elevation (-h <= z <= 0) [m]

    Returns:
    - The value of cosh(k * (z + h)) / cosh(k * h) or an alternative formulation
      for large k * h to prevent floating-point overflow.
    """
    if k * h > 89.4
        # Large k*h: use alternative formulation to avoid overflow
        return exp(k * z) + exp(-k * (z + 2.0 * h))
    else
        # Standard shallow water formulation
        return cosh(k * (z + h)) / cosh(k * h)
    end
end

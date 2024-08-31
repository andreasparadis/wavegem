function one_side_asp(x,t)
# Amplitude spectrum by fast Fourier transform of a timeseries
# x: vector of signal values
# t: time vector

    L = length(x) # Length of the time signal

    T = (t[end]-t[1])/(L-1) # Sampling period of the signal
    fs = 1/T # Sampling frequency

    Nfft = nextpow(2,L)
    y = zeros(Float64,Nfft)
    y[1:L] = x[:]

    H = fft(real(y)) # FFT of the signal
    freq2 = fftfreq(Nfft,fs)  # frequency range
    P2 = abs.(H)/Nfft        # two sided spectrum
    ϕ₂ = angle.(H)

    # single sided spectrum 
    L1 = Int(Nfft/2)+1
    f1 = range(0,fs/2,L1)
    freq1 = collect(f1)
    df = (f1[end]-f1[1])/(length(f1)-1)

    P1 = P2[1:L1]       # single sided spectrum 
    P1[2:end-1] = 2*P1[2:end-1] # f > 0 

    ϕ₁ = ϕ₂[1:L1]       

    return freq1, P1, ϕ₁, df, T, Nfft, H
end

function full_fft_asp!(x,t)
# Amplitude spectrum by fast Fourier transform of a timeseries
# x: vector of signal values
# t: time vector

    L = length(x) # Length of the time signal

    T = (t[end]-t[1])/(L-1) # Sampling period of the signal
    fs = 1/T # Sampling frequency

    # Checks if odd No. of sampling points, if true linearly appends the signals
    t, x, L = fix_odd!(t,x,L,T)

    H = fft(real(x)) # FFT of the signal
    freq2 = fftfreq(L,fs)  # frequency range
    P2 = abs.(H)/L        # two sided spectrum
    ϕ₂ = angle.(H)

    # single sided spectrum 
    L1 = Int(L/2+1)
    f1 = range(0,fs/2,L1)
    freq1 = collect(f1)
    df = (f1[end]-f1[1])/(length(f1)-1)

    P1 = P2[1:L1]       # single sided spectrum 
    P1[2:end-1] = 2*P1[2:end-1] # f > 0 

    ϕ₁ = ϕ₂[1:L1]   

    return H, freq1, P1, ϕ₁, freq2, P2, ϕ₂, x, t, T, L
end

function fix_odd!(t,x,L,T)
    if mod(L,2)==1
        curr_time = Dates.format(now(), "HH:MM:SS")
        print("--------------------------------------------------------\n")
        print("FFT WARNING:   HH:MM:SS = ", curr_time, "               \n")
        print(" - The signal has odd number of elements.               \n")
        print(" - Relevant vectors have been appended linearly,        \n")
        print("   for FFT combatibility.                               \n")
        print(" - Used scheme assumes dt=const.                        \n")
        print("--------------------------------------------------------\n")
        # Append signal with linear interpolation
        push!(t,t[end]+T)
        push!(x,2*x[end]-x[end-1])
        L = L+1
    end

    return t, x, L
end

function phase_shift(arg)
    N=0

    if abs(arg+2*pi) < abs(arg)
        while abs(arg+2*pi*(N+1)) < abs(arg+2*pi*N)
            N=N+1
         end
    elseif abs(arg-2*pi) < abs(arg)
        while abs(arg+2*pi*(N-1)) < abs(arg+2*pi*N)
            N=N-1
        end
    end

    return N
end

function UpDownCross(intime, indata, cs)
    indata[:, 1] .-= mean(indata[:, 1])

    if cs == 0
        iD = findall(diff(sign.(indata[:, 1])) .== 2)  # Up-crossing
    elseif cs == 1
        iD = findall(diff(sign.(indata[:, 1])) .== -2)  # Zero-Down-crossing
    end

    H = Float64[]
    T = Float64[]
    
    for ik in 1:length(iD)-1
        Hk = maximum(indata[iD[ik]:iD[ik+1], 1]) - minimum(indata[iD[ik]:iD[ik+1], 1])
        Tk = maximum(intime[iD[ik]:iD[ik+1], 1]) - minimum(intime[iD[ik]:iD[ik+1], 1])
        push!(H, Hk)
        push!(T, Tk)
    end

    return H, T, iD
end

function decomp_3rd(t0, A0, A1, A2, A3)
    # Obtain spectra for each realization via FFT
    FR1, Â₀₀, ϕ₀₀, _, Tₛ, L = one_side_asp(A0, t0)
    _, Â₀₅, ϕ₀₅, _, _, _ = one_side_asp(A1, t0)
    _, Â₁₀, ϕ₁₀, _, _, _ = one_side_asp(A2, t0)
    _, Â₁₅, ϕ₁₅, _, _, _ = one_side_asp(A3, t0)

    SPEC0 = [FR1 Â₀₀ ϕ₀₀]

    # Spectral amplitudes (not scaled)
    amp00 = Â₀₀*L
    amp10 = Â₀₅*L
    amp05 = Â₁₀*L
    amp15 = Â₁₅*L

    # Spectral phase angles
    arg00 = unwrap(ϕ₀₀)
    arg05 = unwrap(ϕ₀₅)
    arg10 = unwrap(ϕ₁₀)
    arg15 = unwrap(ϕ₁₅)

    # Phase shifting
    sp_peak = findmax(Â₀₀);     iₚ = sp_peak[2];    fₚ = FR1[iₚ]

    N00 = phase_shift(arg00[iₚ])
    arg00 .+= 2 * π * N00
    N05 = phase_shift(arg05[iₚ])
    arg05 .+= 2 * π * N05
    N10 = phase_shift(arg10[iₚ])
    arg10 .+= 2 * π * N10
    N15 = phase_shift(arg15[iₚ])
    arg15 .+= 2 * π * N15

    # Decomposition
    c00 = amp00 .* exp.(1im .* arg00)
    c05 = amp05 .* exp.(1im .* arg05)
    c10 = amp10 .* exp.(1im .* arg10)
    c15 = amp15 .* exp.(1im .* arg15)

    # Amplitudes and phase angles of each component
    ## Linear
    cc1 = ((c00 - c10) - 1im * (c05 - c15)) / 4
    AMP1, ARG1 = abs.(cc1)/L, unwrap(angle.(cc1))
    ## 2nd+
    cc2 = ((c00 + c10) - (c05 + c15)) / 4
    AMP2, ARG2 = abs.(cc2)/L, unwrap(angle.(cc2))
    ## 3rd
    cc3 = ((c00 - c10) + 1im * (c05 - c15)) / 4
    AMP3, ARG3 = abs.(cc3)/L, unwrap(angle.(cc3))
    ## 2nd-
    cc4 = (c00 + c05 + c10 + c15) / 4
    AMP4, ARG4 = abs.(cc4)/L, unwrap(angle.(cc4))

    AMPS = [AMP1 AMP4 AMP2 AMP3]
    ARGS = [ARG1 ARG4 ARG2 ARG3]

    # Obtain time signals of the components via inverse FFT
    L2 = Int(round(L/2))-1;     pad = zeros(Float64,L2)
    lin = ifft([cc1; pad]);     lin = real(lin[1:nₜ])
    plus2 = ifft([cc2; pad]);   plus2 = real(plus2[1:nₜ])
    third = ifft([cc3; pad]);   third = real(third[1:nₜ])
    minus2 = ifft([cc4; pad]);  minus2 = real(minus2[1:nₜ])

    sig_comps = [lin minus2 plus2 third]
    sig_recon = real(lin+plus2+third+minus2)

    return SPEC0, AMPS, ARGS, sig_comps, sig_recon, L, fₚ
end

function low_pass_filter(x,fcut,fsamp,Ntaps)
    responsetype = Lowpass(fcut; fs=fsamp)
    designmethod = FIRWindow(hanning(Ntaps; zerophase=false))
    xₗₚ = filtfilt(digitalfilter(responsetype, designmethod), x)

    return xₗₚ
end
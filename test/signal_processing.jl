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

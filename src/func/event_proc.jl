function ev_identify!(x, t, PeakVal, PeakPos, PeakId, t_range, NP, shift)
    # shift = 0, 1 or 2 == Shift peak at t=0, Start at t=0, No shift
    Lₜ = length(t)
    Tₛ = (t[end]-t[1])/(Lₜ+1)

    PId = sortperm(PeakVal, rev=true)
    PeakVal = PeakVal[PId]
    PeakPos = PeakPos[PId]
    PeakId = PeakId[PId]

    Peaks_len = length(PeakId)

    # Identify and exctract events based on peaks
    rng = round(Int, (t_range/2) / Tₛ + 1)

    ## Move peaks with time range outside the time limits, at the end of the pecking order
    NPout = 0 # No of peaks outside time limits
    for i ∈ 2:Peaks_len-1
        if PeakId[i] < rng+1 || PeakId[i]+rng > Lₜ
            PeakId = [PeakId[1:i-1]; PeakId[i+1:end]; PeakId[i]]
            PeakPos = [PeakPos[1:i-1]; PeakPos[i+1:end]; PeakPos[i]]
            PeakVal = [PeakVal[1:i-1]; PeakVal[i+1:end]; PeakVal[i]]
            NPout += 1 
        end
    end

    ## Isolate events around peaks
    NoPeaks = NP-NPout
    println("$NoPeaks number of qualified events have been identified.")
    
    ∅ = zeros(2*rng+1,NoPeaks)
    tₑᵥ, events = (∅[:,:] for _ = 1:2)

    for i ∈ 1:NoPeaks
        # Time vector of t_range around peak i
        if shift == 0
            tₑᵥ[:, i] = (t[PeakId[i]-rng:PeakId[i]+rng]) .- t[PeakId[i]]
        elseif shift == 1
            tₑᵥ[:, i] = (t[PeakId[i]-rng:PeakId[i]+rng]) .- t[PeakId[i]-rng]
        elseif shift == 2
            tₑᵥ[:, i] = (t[PeakId[i]-rng:PeakId[i]+rng])
        else
            error("Pick a valid shift flag (0, 1 or 2) for the time vectors of events.")
        end

        events[:, i] = x[PeakId[i]-rng:PeakId[i]+rng]
    end

    return tₑᵥ, events, rng, PeakId, PeakPos, PeakVal, NoPeaks
end

function ev_identify_itp!(x, t, PeakVal, PeakPos, PeakId, t_range, NP, shift, res_f)
    # shift = 0, 1 or 2 == Shift peak at t=0, Start at t=0, No shift
    Lₜ = length(t)
    Tₛ = (t[end]-t[1])/(Lₜ+1)

    PId = sortperm(PeakVal, rev=true)
    PeakVal = PeakVal[PId]
    PeakPos = PeakPos[PId]
    PeakId = PeakId[PId]

    Peaks_len = length(PeakId)

    # Identify and exctract events based on peaks
    rng = round(Int, (t_range/2) / Tₛ + 1)

    ## Move peaks with time range outside the time limits, at the end of the pecking order
    NPout = 0 # No of peaks outside time limits
    for i ∈ 2:Peaks_len-1
        if PeakId[i] < rng+1 || PeakId[i]+rng > Lₜ
            PeakId = [PeakId[1:i-1]; PeakId[i+1:end]; PeakId[i]]
            PeakPos = [PeakPos[1:i-1]; PeakPos[i+1:end]; PeakPos[i]]
            PeakVal = [PeakVal[1:i-1]; PeakVal[i+1:end]; PeakVal[i]]
            NPout += 1 
        end
    end

    ## Isolate events around peaks
    NoPeaks = NP-NPout
    println("$NoPeaks number of qualified events have been identified.")
    
    ∅ = zeros(2*rng+1,NoPeaks)
    tₑᵥ, events = (∅[:,:] for _ = 1:2)

    for i ∈ 1:NoPeaks
        # Time vector of t_range around peak i
        if shift == 0
            tₑᵥ[:, i] = (t[PeakId[i]-rng:PeakId[i]+rng]) .- t[PeakId[i]]
        elseif shift == 1
            tₑᵥ[:, i] = (t[PeakId[i]-rng:PeakId[i]+rng]) .- t[PeakId[i]-rng]
        elseif shift == 2
            tₑᵥ[:, i] = (t[PeakId[i]-rng:PeakId[i]+rng])
        else
            error("Pick a valid shift flag (0, 1 or 2) for the time vectors of events.")
        end

        events[:, i] = x[PeakId[i]-rng:PeakId[i]+rng]
    end

    if res_f == 1
        Tₛᵉᵛ = Tₛ/res_f         # Increase resolution by factor res_f
        tₑᵥⁱ = tₑᵥ
        ev_int = events
    else
        # Interpolate events to increase resolution of resulting spectra
        Tₛᵉᵛ = Tₛ/res_f         # Increase resolution by factor res_f
        Lₑᵥ = (2*rng+1) * res_f # Increase length by factor res_f

        ∅ⁱ = zeros(Lₑᵥ, NoPeaks)
        tₑᵥⁱ, ev_int = (∅ⁱ[:,:] for _ = 1:2)

        for i ∈ 1:NoPeaks
            itp = interpolate(tₑᵥ[:, i], events[:, i], BSplineOrder(4))
            tₑᵥⁱ[:, i] = range(tₑᵥ[1],tₑᵥ[end], Lₑᵥ)
            ev_int[:, i] = itp.(tₑᵥⁱ[:, i])
        end
    end

    return tₑᵥⁱ, ev_int, tₑᵥ, events, rng, Tₛᵉᵛ, PeakId, PeakPos, PeakVal, NoPeaks
end

function sect_corr(t_rng, dt, Sects)
    id_rng = round(Int, (t_rng/2) / dt + 1)  # +- in terms of time intervals
    ID0 = div(size(Sects, 1), 2)  # Index of t=0 (middle of vector)

    ## Pulses corresponding to the shortened time range
    shrt = Sects[ID0-id_rng:ID0+id_rng, :]

    # Find the pulse with the highest peak
    M_pls_i = findmax(shrt)[2]
    fmf, inx = M_pls_i[1], M_pls_i[2]

    ## 'inx' is the identifier of the pulse containing the overall max
    A = shrt[:, inx]

    ## Statistical comparison of other pulses against the identified pulse
    ## Exclude pulse A from the set of pulses
    B = hcat(shrt[:, 1:inx-1], shrt[:, inx+1:end])

    cvar = cov(A, B)
    sigA = std(A)
    sigB = std(B, dims=1)

    ## Pearson correlation coefficient
    rho = cvar ./ (sigA * sigB)

    R = A .- B
    Ravg = mean(R, dims=1)
    sigR = std(R, dims=1)

    return rho, Ravg, sigR
end

function instant_freq(x,t)
    # Hilbert transform
    𝓗 = hilbert(x)     
    u = abs.(𝓗)             # envelope
    θ = unwrap(angle.(𝓗))   # phase angle

    # Calculation of instantaneous frequency 
    nₜ = length(t)
    ωᴴ = zeros(Float64, nₜ)     # Instantaneous frequency
    evsep = zeros(Float64, nₜ)  # Event separator line

    for i ∈ 2:nₜ-1
        ωᴴ[i] = (θ[i+1] - θ[i-1]) / (t[i+1]-t[i-1])
        if ωᴴ[i] > 0 
            evsep[i] = maximum(x)
        else
            evsep[i] = 0
        end
    end

    # Statistics of instantaneous frequency
    ω̅ᵢₙₛₜ = mean(ωᴴ)
    idᵤₚ = findall(diff(sign.(ωᴴ)) .== 2) # Up-crossings of instantaneous frequency
    Nᵤₚ = length(idᵤₚ)
    ΔTᵉ = diff(t[idᵤₚ]);          ΔTᵉ = [t[idᵤₚ[1]]; ΔTᵉ]
    tΔT = zeros(Float64,Nᵤₚ);       tΔT[1] = t[idᵤₚ[1]]/2
    for i ∈ 2:Nᵤₚ
        tΔT[i] = (t[idᵤₚ[i-1]] + t[idᵤₚ[i]])/2
    end

    evAmax = zeros(Float64,Nᵤₚ)  
    evAmax[1] = maximum(x[1:idᵤₚ[1]])

    for i ∈ 2:Nᵤₚ
        evAmax[i] = maximum(x[idᵤₚ[i-1]:idᵤₚ[i]])
    end

    return ωᴴ, evsep, ω̅ᵢₙₛₜ, ΔTᵉ, tΔT, idᵤₚ, Nᵤₚ, evAmax
end
function ev_identify(x, t, MinPeakVal, MinPeakDist, t_range, NP, shift, res_f)
    # shift = 0, 1 or 2 == No shift, Start at t=0, Shift peak at t=0
    Lₜ = length(t)
    Tₛ = t[2]-t[1]

    # Find peaks and sort them in descending order
    PeakVal, PeakPos, PeakId = peaks_max_ext(x, t, MinPeakVal, MinPeakDist)

    PId = sortperm(PeakVal, rev=true)
    PeakVal = PeakVal[PId]
    PeakPos = PeakPos[PId]
    PeakId = PeakId[PId]

    Peaks_len = length(PeakId)

    # Identify and exctract events based on peaks
    rng = round(Int, (t_range/2) / Tₛ + 1)

    ## Move peaks with time range outside the time limits, at the end of the pecking order
    NPout = 0 # No of peaks outside time limits
    for i ∈ 1:Peaks_len
        if PeakId[i] < rng+1 || PeakId[i]+rng > Lₜ
            PeakId[NP+1], PeakId[1] = PeakId[1], PeakId[NP+1]
            PeakPos[NP+1], PeakPos[1] = PeakPos[1], PeakPos[NP+1]
            PeakVal[NP+1], PeakVal[1] = PeakVal[1], PeakVal[NP+1]
            NPout += 1 
        end
    end

    ## Isolate events around peaks
    NoPeaks = NP-NPout
    ∅ = zeros(2*rng+1,NoPeaks)
    tₑᵥ, events = (∅[:,:] for _ = 1:2)

    for i ∈ 1:NoPeaks
        # Time vector of t_range around peak i
        if shift == 0
            tₑᵥ[:, i] = (t[PeakId[i]-rng:PeakId[i]+rng])
        elseif shift == 1
            tₑᵥ[:, i] = (t[PeakId[i]-rng:PeakId[i]+rng]) .- t[PeakId[i]-rng]
        elseif shift == 2
            tₑᵥ[:, i] = (t[PeakId[i]-rng:PeakId[i]+rng]) .- t[PeakId[i]]
        else
            error("Pick a valid shift flag (0, 1 or 2) for the time vectors of events.")
        end

        events[:, i] = x[PeakId[i]-rng:PeakId[i]+rng]
    end

    println("$NoPeaks number of qualified events have been identified.")

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

    return tₑᵥⁱ, ev_int, tₑᵥ, events, rng, Tₛᵉᵛ, PeakId, PeakPos, PeakVal
end
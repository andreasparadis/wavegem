# Peak detection algorithm based on the 2nd derivative of the input
# signal, using a 1st order symmetric finite difference scheme

function peaks(x,t)
    # Initialize variables
    m = length(t)
    𝚶⃗ₘ = zeros(Float64,m);    𝚶⃗₁ = zeros(Float64,1)
    ẋ, ẍ = 𝚶⃗ₘ[:], 𝚶⃗ₘ[:]
    t⁺, t¯, x⁺, x¯ = (𝚶⃗₁[:] for _ = 1:4)
    i⁺ = zeros(Int64,1);    i¯ = zeros(Int64,1)

    # 1st differentiation
    for i ∈ 2:m-1
        ẋ[i] = (x[i+1] - x[i-1]) / (t[i+1]-t[i-1])
    end

    # Linear interpolation for first element
    ẋ[1] = 2*ẋ[2]-ẋ[3]

    # 2nd differentiation and storing of maxima and minima
    for i ∈ 2:m-1
        ẍ[i] = (ẋ[i+1] - ẋ[i-1]) / (t[i+1]-t[i-1])
        if sign(ẋ[i]) != sign(ẋ[i-1])
            if ẍ[i-1] < 0
                push!(i⁺,i-1)
                push!(t⁺,t[i-1])
                push!(x⁺,x[i-1])
            else
                push!(i¯,i-1)
                push!(t¯,t[i-1])
                push!(x¯,x[i-1])
            end
        end
    end

    # Linear interpolation of last element
    ẍ[end] = 2*ẍ[end-1]-ẍ[end-2]

    i⁺ = i⁺[2:end];   i¯ = i¯[2:end]
    t⁺ = t⁺[2:end];   t¯ = t¯[2:end]
    x⁺ = x⁺[2:end];   x¯ = x¯[2:end]
    return t⁺, t¯, x⁺, x¯, ẋ, ẍ, i⁺, i¯
end

function peaks_extend(x,t, MinPeakVal, MinPeakDist)
    # Initialize variables
    m = length(t)
    𝚶⃗ₘ = zeros(Float64,m);    
    Ø = Array{Float64}(undef,0)
    ẋ, ẍ = 𝚶⃗ₘ[:], 𝚶⃗ₘ[:]
    PeakVal, PeakPos = Ø[:], Ø[:]
    PeakId = Array{Int64}(undef,0)   

    # 1st differentiation
    for i ∈ 2:m-1
        ẋ[i] = (x[i+1] - x[i-1]) / (t[i+1]-t[i-1])
    end

    # Linear interpolation for first element
    ẋ[1] = 2*ẋ[2]-ẋ[3]

    # 2nd differentiation and storing of maxima and minima
    for i ∈ 2:m-1
        ẍ[i] = (ẋ[i+1] - ẋ[i-1]) / (t[i+1]-t[i-1])
        if sign(ẋ[i]) != sign(ẋ[i-1])
            push!(PeakId,i-1)
            push!(PeakPos,t[i-1])
            push!(PeakVal,x[i-1])
        end
    end

    # Linear interpolation of last element
    ẍ[end] = 2*ẍ[end-1]-ẍ[end-2]

    # Apply the MinPeakVal, MinPeakDist conditions and create the vectors corresponding to these extrema
    # MinPeakVal condition
    ids = findall(f -> abs(f)>MinPeakVal, PeakVal) 
    PId = PeakId[ids]
    PPos = PeakPos[ids]
    PVal = PeakVal[ids]

    n = length(ids)
    tₑ, xₑ = Ø, Ø
    ID = Array{Int64}(undef,0)
    tₑ, xₑ, ID = [PPos[1]], [PVal[1]], [PId[1]]

    # MinPeakDist condition
    for i ∈ 2:n
        if PPos[i]-PPos[i-1] > MinPeakDist
            push!(ID, PId[i])
            push!(tₑ, PPos[i])
            push!(xₑ, PVal[i])
        end
    end

    return xₑ, tₑ, ID, ẋ, ẍ
end

function maxima(x,t)
    # Initialize variables
    m = length(t)
    𝚶⃗ₘ = zeros(Float64,m);    𝚶⃗₁ = zeros(Float64,1)
    ẋ, ẍ = 𝚶⃗ₘ[:], 𝚶⃗ₘ[:]
    t⁺ = 𝚶⃗₁[:];  x⁺ = 𝚶⃗₁[:];    i⁺ = zeros(Int64,1)

    # 1st differentiation
    for i ∈ 2:m-1
        ẋ[i] = (x[i+1] - x[i-1]) / (t[i+1]-t[i-1])
    end

    # Linear interpolation of last element
    ẋ[1] = 2*ẋ[2]-ẋ[3]

    # 2nd differentiation and storing of maxima and minima
    for i ∈ 2:m-1
        ẍ[i] = (ẋ[i+1] - ẋ[i-1]) / (t[i+1]-t[i-1])
        if sign(ẋ[i]) != sign(ẋ[i-1])
            if ẍ[i-1] < 0
                push!(i⁺,i-1)
                push!(t⁺,t[i-1])
                push!(x⁺,x[i-1])
            end
        end
    end

    # Linear interpolation of last element
    ẍ[end] = 2*ẍ[end-1]-ẍ[end-2]

    i⁺ = i⁺[2:end]
    t⁺ = t⁺[2:end]
    x⁺ = x⁺[2:end]
    return t⁺, x⁺, ẋ, ẍ, i⁺
end

function peaks_max_ext(x,t, MinPeakVal, MinPeakDist, sort)
    # x: signal ,t: time, MinPeakVal, MinPeakDist, sort: true/false - sort peaks by amplitude
    # Initialize variables
    m = length(t)
    𝚶⃗ₘ = zeros(Float64,m);    
    Ø = Array{Float64}(undef,0)
    ẋ, ẍ = 𝚶⃗ₘ[:], 𝚶⃗ₘ[:]
    PeakVal, PeakPos = Ø[:], Ø[:]
    PeakId = Array{Int64}(undef,0)   

    # 1st differentiation
    for i ∈ 2:m-1
        ẋ[i] = (x[i+1] - x[i-1]) / (t[i+1]-t[i-1])
    end

    # Linear interpolation for first element
    ẋ[1] = 2*ẋ[2]-ẋ[3]

    # 2nd differentiation and storing of maxima and minima
    for i ∈ 2:m-1
        ẍ[i] = (ẋ[i+1] - ẋ[i-1]) / (t[i+1]-t[i-1])
        if sign(ẋ[i]) != sign(ẋ[i-1]) && ẍ[i-1] < 0
            push!(PeakId,i-1)
            push!(PeakPos,t[i-1])
            push!(PeakVal,x[i-1])
        end
    end

    # Linear interpolation of last element
    ẍ[end] = 2*ẍ[end-1]-ẍ[end-2]

    # Apply the MinPeakVal, MinPeakDist conditions and create the vectors corresponding to these extrema
    # MinPeakVal condition
    ids = findall(f -> f>MinPeakVal, PeakVal) 
    PId = PeakId[ids]
    PPos = PeakPos[ids]
    PVal = PeakVal[ids]

    n = length(ids)
    tₑ, xₑ = Ø, Ø
    ID = Array{Int64}(undef,0)

    if sort
        # Sort in descending order based on amplitude
        for i in 1:n
            for j in 1:n-i
                if PVal[j] < PVal[j+1]
                    PVal[j], PVal[j+1] = PVal[j+1], PVal[j]
                    PPos[j], PPos[j+1] = PPos[j+1], PPos[j]
                    PId[j], PId[j+1] = PId[j+1], PId[j]
                end
            end
        end
    end

    tₑ, xₑ, ID = [PPos[1]], [PVal[1]], [PId[1]]

    # MinPeakDist condition
    # cnt = 1
    for i ∈ 2:n
        PeakDist = abs.(tₑ .- PPos[i])
        if minimum(PeakDist) > MinPeakDist
            push!(ID, PId[i])
            push!(tₑ, PPos[i])
            push!(xₑ, PVal[i])
            # cnt = i
        end
    end

    return xₑ, tₑ, ID, ẋ, ẍ
end
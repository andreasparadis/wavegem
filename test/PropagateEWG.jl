# Propagation of individual EWGs
TOT = zeros(ComplexF64, M)  # Total surface elevation
Nₛ = Int(nextpow(2,M)/2)+1
SP_TOT = zeros(Float64, Nₛ) # Total spectrum of propagated Gaussian EWGs

pltEWG = plot(palette=:darkrainbow)
pltSPEC = plot(palette=:darkrainbow)
pltCONV = plot(palette=:darkrainbow)

for n ∈ 1:lenT
    # ξ = -β .+ [π/2; π; π; π/2] # .+ zeros(Float64, lenT)
    # ξ = ωᵢ[n] * t .+ β
     
    # Gaussian EWG envelope
    gₙ = elem_wg(t̅*Tₚⱼ, A̅ₒ*Hₛⱼ, t̅ᶜ*Tₚⱼ, T̅ₒ*Tₚⱼ, n)  # Non-dimensional
    # gₙ = elem_wg(t, Aₒ, tᶜₒ, Tₒ, n)

    # Propagated EWG
    WGT,_ = wave_group(1,T̅ₒ[n]*Tₚⱼ,1,2π/(Ω*ωₚⱼ),d,t̅*Tₚⱼ,0,β̃[n])
    # # WGT,_ = wave_group(1,T̅ₒ[n]*Tₚⱼ,1,2π/(ωₚ),d,t̅*Tₚⱼ,0,0)  
    # # WGT,_ = stokes_sum(1,T̅ₒ[n]*Tₚⱼ,1,2π/(Ω*ωₚⱼ),d,t̅*Tₚⱼ,0,β̃[n])
    # ηₙ = gₙ .* WGT 
    # # ηₙ = gₙ .* exp.(-1im * Ω*ωₚⱼ * t̅*Tₚⱼ) #.* exp.(1im * β)
    # # ηₙ = gₙ .* exp.(-1im * ωₒ[n] * t) .* exp.(1im * β)
    # FR, MAG, ang, df,_ = one_side_asp(real(ηₙ),t̅*Tₚⱼ) 

    # Propagated EWG - Focused WG
    ηⁱₙ = gₙ .* exp.(-1im * Ω*ωₚⱼ * t̅*Tₚⱼ) 
    # ηⁱₙ = gₙ .* exp.(-1im * ωₚⱼ * t̅*Tₚⱼ) 
    # ηⁱₙ = gₙ .* exp.(-1im * ωₒ[n] * t̅*Tₚⱼ)
    FR, MAG, ang, df,_ = one_side_asp(real(ηⁱₙ),t̅*Tₚⱼ)
    ηₙ = fcsd_wg(FR, MAG, t̅*Tₚⱼ, Aₒ[n], t̅ᶜ[n]*Tₚⱼ, β̃[n])

    # Focused WG - JONSWAP
    # FR,MAG = spectrum(Aₒ[n],2π/ω̃,γ,Tcut,128,Nₛ) # Generate JONSWAP spectrum
    # ηₙ = fcsd_wg(FR, MAG, t̅*Tₚⱼ, Aₒ[n], t̅ᶜ[n]*Tₚⱼ, β̃[n])       

    TOT[:] = TOT[:] .+ ηₙ
    SP_TOT[:] = SP_TOT[:] .+ MAG

    plot!(pltEWG, t̅*Tₚⱼ, real(ηₙ), xlab = "t [s]", ylab = "A [m]", lab = "g$(n)", lw=2, legend=:outerbottom, legendcolumns=lenT)
    plot!(pltSPEC, FR, MAG, xlab = L"f [Hz]", ylab = L"A [m^2 s]", lab = "g$(n)", lw=2, legend=:outerbottom, legendcolumns=lenT)
    plot!(pltCONV, Tₒₚₜ[2][:,n], xlab = "iteration",  ylab = L"T_o^i [s]", xscale=:log10, lab = "To[$(n)]", lw=2, legend=:outerbottom, legendcolumns=lenT)

    for i ∈ 1:M
        if abs(ηₙ[i]) < 1e-6
            ηₙ[i] = 0
        end
    end                             

    fid_EWG = joinpath(fid_res,"EWG_$n") # File name
    open(fid_EWG, "w")
    # writedlm(fid_EWG, [round.(t*1e6)./1e6 round.(gₙ*1e6)./1e6 round.(real(ηₙ*1e6))./1e6 round.(angle.(ηₙ)*1e6)./1e6], '\t')
    writedlm(fid_EWG, [round.(t*1e6)./1e6 round.(gₙ*1e6)./1e6 round.(WGT*1e6)./1e6 round.(ηₙ*1e6)./1e6], '\t')
end
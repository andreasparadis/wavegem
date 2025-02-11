function insta(Px,Pz,Qx,Qz)
    N = length(Px)
    xC, zC = zeros(Float64,N), zeros(Float64,N)

    for i ∈ 1:N-1
        P1 = [Px[i];Pz[i]]
        Q1 = [Qx[i];Qz[i]]
        P2 = [Px[i+1];Pz[i+1]]
        Q2 = [Qx[i+1];Qz[i+1]]

        # P1 = [189.03668300000024;-281.1793970000003]
        # Q1 = [123.93971000000012;-280.8299280000003]
        # P2 = [189.0229464376222;-281.1818961242428]
        # Q2 = [123.92867661384555;-280.82621261044665]

        # P1 = [X[3000,1];Z[3000,1]]
        # Q1 = [X[3000,4];Z[3000,4]]
        # P2 = [X[3001,1];Z[3001,1]]
        # Q2 = [X[3001,4];Z[3001,4]]

        l1 = (P2[2]-P1[2])/(P2[1]-P1[1])
        l2 = (Q2[2]-Q1[2])/(Q2[1]-Q1[1])

        b1 = P1[2]-l1*P1[1]
        b2 = Q1[2]-l2*Q1[1]

        r1 = sqrt((P1[1]-P2[1])^2 + (P1[2]-P2[2])^2) / 2
        r2 = sqrt((Q1[1]-Q2[1])^2 + (Q1[2]-Q2[2])^2) / 2

        A1 = 1+l1^2
        B1 = -2*(P1[1]-l1*(b1-P1[2]))
        C1 = P1[1]^2 + (b1-P1[2])^2-r1^2
        D1 = B1^2-4*A1*C1

        xP12 = (-B1+sqrt(D1))/(2*A1)
        if sign(xP12) != sign(Px[i])
            xP12 = (-B1-sqrt(D1))/(2*A1)
        end
        zP12 = l1*xP12 + b1

        A2 = 1+l2^2
        B2 = -2*(Q1[1]-l2*(b2-Q1[2]))
        C2 = Q1[1]^2 + (b2-Q1[2])^2-r2^2
        D2 = B2^2-4*A2*C2

        xQ12 = (-B2+sqrt(D2))/(2*A2)
        if sign(xQ12) != sign(Qx[i])
            xQ12 = (-B2-sqrt(D2))/(2*A2)
        end
        zQ12 = l2*xQ12 + b2

        ba = zP12+l1*xP12
        bb = zQ12+l2*xQ12

        XC = (bb-ba)/(l2-l1)
        ZC = -l1*XC+ba

        # xC[i] = P1[1]-XC
        # zC[i] = P1[2]-ZC
        xC[i] = XC
        zC[i] = ZC

        # println("xP1-xC = ", P1[1]-xC)
        # println("zP1-zC = ", P1[2]-zC)
    end

    plt_move = plot(xlab="X [mm]", ylab="Z [mm]")
    for i ∈ 1:100:N
        plot!([Px[i];Qx[i]], [Pz[i];Qz[i]], lab="t=$(round(i/120))")
    end

    display(plt_move)
    return xC, zC
end
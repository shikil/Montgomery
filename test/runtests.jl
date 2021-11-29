using Test
using Montgomery

@testset "SIKE for beginners" begin
    # Public parameters
    k, i = FiniteImaginaryQuadraticField(431)
    Ea0 = MontgomeryCurve(329i + 423)
    Pa = Ea0(100i + 248, 304i + 199)
    Qa = Ea0(426i + 394, 51i + 79)
    Pb = Ea0(358i + 275, 410i + 104)
    Qb = Ea0(20i + 185, 281i + 239)
    @test 2Pa == Pa * 2
    @test (100i + 248, 304i + 199) in Ea0
    @test issuper_singular(Ea0)

    # Alice's public key generation
    ka = 11
    Sa = Pa + ka * Qa
    @test Sa.coordx == 271i + 79 && Sa.coordy == 153i + 430
    jvalue = (107, 344i + 190, 350i + 65, 222i + 118)
    PKa = ()
    for (index, value) in enumerate(3:-1:0)
        Ra = 2^value * Sa
        I = Isogeny(Ra, 2)
        Ea = codomain(I)
        @test j_invariant(Ea) == jvalue[index]
        Pb = I(Pb)
        Qb = I(Qb)
        if value > 0
            Sa = I(Sa)
            @test order(Sa) == 2^value
        else
            PKa = (Ea, Pb, Qb)
            @test Pb.coordx == 142i + 183 && Qb.coordx == 220i + 314
        end
    end

    # Bob's public key generation
    kb = 2
    Pb = Ea0(358i + 275, 410i + 104)
    Qb = Ea0(20i + 185, 281i + 239)
    Sb = Pb + kb * Qb
    @test order(Sb) == 27
    jvalue = (106i + 379, 325i + 379, 344i + 190)
    PKb = ()
    for (index, value) in enumerate(2:-1:0)
        Rb = 3^value * Sb
        I = Isogeny(Rb, 3)
        Ea = codomain(I)
        @test j_invariant(Ea) == jvalue[index]
        Pa = I(Pa)
        Qa = I(Qa)
        if value > 0
            Sb = I(Sb)
            @test order(Sb) == 3^value
        else
            PKb = (Ea, Pa, Qa)
            @test Pa.coordx == 187i + 226 && Qa.coordx == 325i + 415
        end
    end

    # Alice's shared secret computation
    # Hard code point value because of general addition and multiplication depends on y-coordinate
    Ea0 = PKb[1]
    Pa = Ea0(187i + 226, 43i + 360)
    Qa = Ea0(325i + 415, 322i + 254)
    Sa = Pa + ka * Qa
    @test Sa.coordx == 125i + 357
    jvalue = (364i + 304, 67, 242, 234)
    E1 = Ea0
    for (index, value) in enumerate(3:-1:0)
        Ra = 2^value * Sa
        I = Isogeny(Ra, 2)
        Ea = codomain(I)
        @test degree(I) == 2
        @test j_invariant(Ea) == jvalue[index]
        if value > 0
            Sa = I(Sa)
        else
            E1 = Ea
        end
    end

    # Bob's shared secret computation
    Ea0 = PKa[1]
    Pb = Ea0(142i + 183, 119i + 360)
    Qb = Ea0(220i + 314, 289i + 10)
    Sb = Pb + kb * Qb
    @test Sb.coordx == 393i + 124
    jvalue = (299i + 315, 61, 234)
    E2 = Ea0
    for (index, value) in enumerate(2:-1:0)
        Rb = 3^value * Sb
        I = Isogeny(Rb, 3)
        Ea = codomain(I)
        @test degree(I) == 3
        @test j_invariant(Ea) == jvalue[index]
        if value > 0 
            Sb = I(Sb)
        else
            E2 = Ea
        end
    end

    @test isisogenous(E1, E2)
end

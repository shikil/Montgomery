using Test
using Montgomery

@testset "SIKE for beginners" begin
    # Public parameters
    k, i = FiniteField(431)
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
    Ra = 8Sa
    Ea1, ϕ = isogeny2(Ra)
    @test j_invariant(Ea1) == 107
    Pb = Ea1(ϕ(Pb.coordx))
    Qb = Ea1(ϕ(Qb.coordx))
    Sa = Ea1(ϕ(Sa.coordx))
    @test order(Sa) == 8

    Ra = 4Sa
    Ea2, ϕ = isogeny2(Ra)
    @test j_invariant(Ea2) == 344i + 190
    Pb = Ea2(ϕ(Pb.coordx))
    Qb = Ea2(ϕ(Qb.coordx))
    Sa = Ea2(ϕ(Sa.coordx))
    @test order(Sa) == 4

    Ra = 2Sa
    Ea3, ϕ = isogeny2(Ra)
    @test j_invariant(Ea3) == 350i + 65
    Pb = Ea3(ϕ(Pb.coordx))
    Qb = Ea3(ϕ(Qb.coordx))
    Sa = Ea3(ϕ(Sa.coordx))
    @test order(Sa) == 2

    Ea4, ϕ = isogeny2(Sa)
    @test j_invariant(Ea4) == 222i + 118
    Pb = Ea4(ϕ(Pb.coordx))
    Qb = Ea4(ϕ(Qb.coordx))
    PKa = (Ea4, Pb, Qb)
    @test Pb.coordx == 142i + 183 && Qb.coordx == 220i + 314

    # Bob's public key generation
    kb = 2
    Pb = Ea0(358i + 275, 410i + 104)
    Qb = Ea0(20i + 185, 281i + 239)
    Sb = Pb + kb * Qb
    @test order(Sb) == 27

    Rb = 9Sb
    Ea1, ϕ = isogeny3(Rb)
    @test j_invariant(Ea1) == 106i + 379
    Pa = Ea1(ϕ(Pa.coordx))
    Qa = Ea1(ϕ(Qa.coordx))
    Sb = Ea1(ϕ(Sb.coordx))
    @test order(Sb) == 9

    Rb = 3Sb
    Ea2, ϕ = isogeny3(Rb)
    @test j_invariant(Ea2) == 325i + 379
    Pa = Ea2(ϕ(Pa.coordx))
    Qa = Ea2(ϕ(Qa.coordx))
    Sb = Ea2(ϕ(Sb.coordx))
    @test order(Sb) == 3

    Ea3, ϕ = isogeny3(Sb)
    @test j_invariant(Ea3) == 344i + 190
    Pa = Ea3(ϕ(Pa.coordx))
    Qa = Ea3(ϕ(Qa.coordx))
    PKb = (Ea3, Pa, Qa)
    @test Pa.coordx == 187i + 226 && Qa.coordx == 325i + 415

    # Alice's shared secret computation
    # Hard code point value because of addition and multiplication depends on y-coordinate
    Ea0 = PKb[1]
    Pa = Ea0(187i + 226, 43i + 360)
    Qa = Ea0(325i + 415, 322i + 254)
    Sa = Pa + ka * Qa
    @test Sa.coordx == 125i + 357
    jvalue = (364i + 304, 67, 242, 234)
    for (index, value) in enumerate(3:-1:0)
        Ra = 2^value * Sa
        Ea, ϕ = isogeny2(Ra)
        @test j_invariant(Ea) == jvalue[index]
        if value > 0
            Sa = Ea(ϕ(Sa.coordx))
        end
    end

    # Bob's shared secret computation
    Ea0 = PKa[1]
    Pb = Ea0(142i + 183, 119i + 360)
    Qb = Ea0(220i + 314, 289i + 10)
    Sb = Pb + kb * Qb
    @test Sb.coordx == 393i + 124
    jvalue = (299i + 315, 61, 234)
    for (index, value) in enumerate(2:-1:0)
        Rb = 3^value * Sb
        Ea, ϕ = isogeny3(Rb)
        @test j_invariant(Ea) == jvalue[index]
        if value > 0 
            Sb = Ea(ϕ(Sb.coordx))
        end
    end
end

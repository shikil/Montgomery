using Test
using Montgomery

@testset "SIKE for beginners" begin
    k, i = FiniteField(431)
    Ea0 = MontgomeryCurve(329i + 423)
    Pa = Ea0(100i + 248, 304i + 199)
    Qa = Ea0(426i + 394, 51i + 79)
    Pb = Ea0(358i + 275, 410i + 104)
    Qb = Ea0(20i + 185, 281i + 239)
    @test 2Pa == Pa * 2
    @test (100i + 248, 304i + 199) in Ea0

    ka = 11
    Sa = Pa + ka * Qa
    @test Sa.coordx == 271i + 79 && Sa.coordy == 153i + 430
    Ra = 8Sa
    Ea1, ϕ = isogeny(Ra)
    @test j_invariant(Ea1) == 107
    @test ϕ(Pb.coordx) == 118i + 85 && ϕ(Qb.coordx) == 62i + 124 && ϕ(Sa.coordx) == 36i + 111
    @test Ea1(118i + 85) in [Ea1(118i + 85, 274i + 150), -Ea1(118i + 85, 274i + 150)]
end

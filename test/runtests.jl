using Test
using Hecke
import Montgomery: MontgomeryCurve

T, t = ResidueRing(ZZ, 431)["t"]
k, i = FiniteField(t^2 + 1, "i")
Ea0 = MontgomeryCurve(329i+423)

@test 2Ea0([100i+248, 304i+199]) == Ea0([100i+248, 304i+199]) * 2
@test Ea0(100i+248) in [Ea0([100i+248, 304i+199]), -Ea0([100i+248, 304i+199])]

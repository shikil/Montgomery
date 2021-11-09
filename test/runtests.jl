using Test
using Montgomery

k, i = FiniteField(431)
Ea0 = MontgomeryCurve(329i+423)

@test 2Ea0([100i+248, 304i+199]) == Ea0([100i+248, 304i+199]) * 2
@test Ea0(100i+248) in [Ea0([100i+248, 304i+199]), -Ea0([100i+248, 304i+199])]

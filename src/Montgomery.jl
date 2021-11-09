module Montgomery

import Base.:*
import Hecke: FieldElem, EllipticCurve, EllCrv, EllCrvPt, ResidueRing, ZZ, j_invariant, parent, root
import Nemo: FlintFiniteField
export FiniteField, MontgomeryCurve, j_invariant

Base.:*(P::EllCrvPt, n::Int) = n * P

function (E::EllCrv{T})(a::T) where T
	E([a, root(a^3 + E.coeff[2] * a^2 + E.coeff[4] * a + E.coeff[5], 2)])
end

function (E::EllCrv{T})(x::T, y::T) where T
	E([x, y])
end

function FiniteField(p; gen="i")
	_, t = ResidueRing(ZZ, p)["t"]
	FlintFiniteField(t^2+1, gen)
end

function MontgomeryCurve(a::T) where T<:FieldElem
    aid = zero(parent(a))
    EllipticCurve([aid, a, aid, one(a), aid])
end

end # module

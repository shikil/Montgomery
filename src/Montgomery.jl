module Montgomery

import Hecke: EllipticCurve, EllCrv, EllCrvPt, ResidueRing, ZZ, ison_curve, j_invariant, order, parent, roots
import Nemo: FlintFiniteField

export FiniteField, MontgomeryCurve, isogeny, j_invariant

Base.:*(p::EllCrvPt, n::Int) = n * p

function Base.in(coords::Tuple{T, T}, e::EllCrv{T}) where T
	ison_curve(e, [coords[1], coords[2]])
end

function (e::EllCrv{T})(a::T) where T
	_, y = parent(a)["y"]
	e(a, roots(y^2 - (a^3 + e.coeff[2] * a^2 + e.coeff[4] * a + e.coeff[5]))[1])
end

function (e::EllCrv{T})(x::T, y::T) where T
	e([x, y])
end

function FiniteField(p; gen="i")
	_, t = ResidueRing(ZZ, p)["t"]
	FlintFiniteField(t^2+1, gen)
end

function MontgomeryCurve(a::T) where T
    aid = zero(parent(a))
    EllipticCurve([aid, a, aid, one(a), aid])
end

function isogeny(p::EllCrvPt)
	a = 2 * (1 - 2p.coordx^2)
	MontgomeryCurve(a), x -> x * (p.coordx * x - 1)//(x - p.coordx)
end

end # module

module Montgomery
import Base.:*
import Hecke: FieldElem, EllipticCurve, EllCrv, EllCrvPt, parent, root
export MontgomeryCurve

Base.:*(P::EllCrvPt, n::Int64) = n * P

function (E::EllCrv{T})(a::T) where T
	E([a, root(a^3 + E.coeff[2] * a^2 + E.coeff[4] * a + E.coeff[5], 2)])
end

function MontgomeryCurve(a::T) where T<:FieldElem
    aid = zero(parent(a))
    EllipticCurve([aid, a, aid, one(a), aid])
end

end # module

module Montgomery

import AbstractAlgebra: FinFieldElem, elem_type
import Hecke: EllipticCurve, EllCrv, EllCrvPt, ResidueRing, ZZ, base_field, characteristic, ison_curve, j_invariant, order, parent, roots
import Nemo: FlintFiniteField

export FiniteImaginaryQuadraticField, MontgomeryCurve, isisogenous, isogeny2, isogeny3, issuper_singular, j_invariant, order

elem_type(::Type{EllCrv{T}}) where T = EllCrvPt{T}

Base.:*(p::EllCrvPt, n::Int) = n * p

Φ(x,y) = -x^2 * y^2 + x^3 + 1488x^2 * y + 1488x * y^2 + y^3 - 162000x^2 + 40773375x * y - 162000y^2 + 8748000000x + 8748000000y - 157464000000000

function Base.in(coords::Tuple{T,T}, e::EllCrv{T}) where T
    ison_curve(e, [coords[1], coords[2]])
end

function (e::EllCrv{T})(a::T) where T
    _, y = parent(a)["y"]
    e(a, roots(y^2 - (a^3 + e.coeff[2] * a^2 + e.coeff[4] * a + e.coeff[5]))[1])
end

function (e::EllCrv{T})(x::T, y::T) where T
    e([x, y])
end

function FiniteImaginaryQuadraticField(p; gen="i")
    _, t = ResidueRing(ZZ, p)["t"]
    FlintFiniteField(t^2 + 1, gen)
end

MontgomeryCurve(a) = EllipticCurve(parent(a).([0, a, 0, 1, 0]))

function isisogenous(e1::EllCrv{T}, e2::EllCrv{T}) where T <: FinFieldElem
    j_invariant(e1) == j_invariant(e2)
end

function isogeny2(p::EllCrvPt)
    a = 2 * (1 - 2p.coordx^2)
    MontgomeryCurve(a), x -> x * (p.coordx * x - 1) // (x - p.coordx)
end

function isogeny3(p::EllCrvPt)
    a = (parent(p).coeff[2] * p.coordx - 6p.coordx^2 + 6) * p.coordx
    MontgomeryCurve(a), x -> x * (p.coordx * x - 1)^2 // (x - p.coordx)^2
end

function issuper_singular(e::EllCrv{T}) where T <: FinFieldElem
    order(e) % characteristic(base_field(e)) == 1
end

end # module

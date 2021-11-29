module Montgomery

import AbstractAlgebra: FinFieldElem, Map, SetMap, codomain, elem_type
import Hecke: EllipticCurve, EllCrv, EllCrvPt, ResidueRing, ZZ, base_field, characteristic, ison_curve, j_invariant, order, parent, roots
import Nemo: FlintFiniteField

export FiniteImaginaryQuadraticField, Isogeny, MontgomeryCurve, codomain, degree, isisogenous, issuper_singular, j_invariant, order

struct Isogeny <: Map{EllCrv, EllCrv, SetMap, Isogeny}
    domain::EllCrv
    codomain::EllCrv
    degree::Int
    _ϕ

    function Isogeny(p::EllCrvPt{T}, degree::Int) where T <: FinFieldElem
        if degree == 2
            codomain = MontgomeryCurve(2(1 - 2p.coordx^2))
            ϕ = x -> x * (p.coordx * x - 1) // (x - p.coordx)
        elseif degree == 3
            codomain = MontgomeryCurve((parent(p).coeff[2] * p.coordx - 6p.coordx^2 + 6) * p.coordx)
            ϕ = x -> x * (p.coordx * x - 1)^2 // (x - p.coordx)^2
        end
        new(parent(p), codomain, degree, ϕ)
    end
end

(i::Isogeny)(x::Int) = i._ϕ(x)
(i::Isogeny)(e::EllCrvPt{T}) where T <: FinFieldElem = codomain(i)(i._ϕ(e.coordx))
degree(f::Isogeny) = getfield(f, :degree)

Base.:*(p::EllCrvPt, n::Int) = n * p
Base.in(coords::Tuple{T,T}, e::EllCrv{T}) where T = ison_curve(e, [coords[1], coords[2]])

elem_type(::Type{EllCrv{T}}) where T = EllCrvPt{T}

Φ(x,y) = -x^2 * y^2 + x^3 + 1488x^2 * y + 1488x * y^2 + y^3 - 162000x^2 + 40773375x * y - 162000y^2 + 8748000000x + 8748000000y - 157464000000000

function (e::EllCrv{T})(a::T) where T
    _, y = parent(a)["y"]
    e(a, roots(y^2 - (a^3 + e.coeff[2] * a^2 + e.coeff[4] * a + e.coeff[5]))[1])
end

(e::EllCrv{T})(x::T, y::T) where T = e([x, y])

function FiniteImaginaryQuadraticField(p; gen="i")
    _, t = ResidueRing(ZZ, p)["t"]
    FlintFiniteField(t^2 + 1, gen)
end

MontgomeryCurve(a) = EllipticCurve(parent(a).([0, a, 0, 1, 0]))

isisogenous(e1::EllCrv{T}, e2::EllCrv{T}) where T <: FinFieldElem = j_invariant(e1) == j_invariant(e2)

issuper_singular(e::EllCrv{T}) where T <: FinFieldElem = order(e) % characteristic(base_field(e)) == 1

end # module

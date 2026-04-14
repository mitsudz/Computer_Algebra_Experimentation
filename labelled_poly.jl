#=
# Fast Labelled Polynomials for F5B
# Uses bit-packed GrLex representation for both signatures and polynomials.
=#

using DataStructures

# Note: Assumes FastPoly, GrLexMonomial, etc., are defined in your environment

# TODO - change index to a UInt16 with the assumption that we don't get 60k+ initial polynomials, also better name for index
# TODO - Move signature and index to their own struct called signature
struct LabelledPolynomial{C}
    index::Int
    signature::GrLexMonomial
    poly::FastPoly{C}
end

"""
Signature comparison logic for F5.
In F5, (S1, i) > (S2, j) if:
1. i < j  (Lower index is 'older' and thus greater in some F5 variations)
2. i == j and S1 > S2 (Standard monomial order comparison)
"""
@inline function is_sig_greater(LP1::LabelledPolynomial, LP2::LabelledPolynomial)
    LP1.index == LP2.index && return LP1.signature > LP2.signature
    return LP1.index < LP2.index
end

# --- Getters ---

@inline leading_monomial(F::LabelledPolynomial) = leading_monomial(F.poly)
@inline leading_term(F::LabelledPolynomial) = leading_term(F.poly)
@inline leading_coefficient(F::LabelledPolynomial) = leading_coefficient(F.poly)

# --- Arithmetic ---

"""
Multiplication of a LabelledPolynomial by a Term or Monomial.
Maintains the F5 property: sig(t * F) = mono(t) * sig(F).
"""
@inline function Base.:*(LP::LabelledPolynomial{C}, t::Term{C}) where C
    return LabelledPolynomial(
        LP.index,
        t.mono * LP.signature,
        t * LP.poly
    )
end
@inline Base.:*(a::Term{C}, b::LabelledPolynomial{C}) where C = b * a

@inline function Base.:*(LP::LabelledPolynomial{C}, m::GrLexMonomial) where C
    return LabelledPolynomial(
        LP.index,
        m * LP.signature,
        m * LP.poly
    )
end
@inline Base.:*(a::GrLexMonomial, b::LabelledPolynomial{C}) where C = b * a

@inline function Base.:*(LP::LabelledPolynomial{C}, c::Number) where C
    # Multiplying by a constant does not change the signature monomial
    return LabelledPolynomial(LP.index, LP.signature, c * LP.poly) 
end
@inline Base.:*(a::Number, b::LabelledPolynomial) = b * a

"""
Addition of two LabelledPolynomials.
The signature of the sum is the maximum of the two input signatures.

Precondition: LP1.signature != LP2.signature (otherwise infinite recursion)
"""
@inline function Base.:+(LP1::LabelledPolynomial{C}, LP2::LabelledPolynomial{C}) where C
    is_sig_greater(LP1, LP2) && return LabelledPolynomial(LP1.index, LP1.signature, LP1.poly + LP2.poly)
    return LabelledPolynomial(LP2.index, LP2.signature, LP1.poly + LP2.poly)
end

@inline function Base.:-(LP1::LabelledPolynomial{C}, LP2::LabelledPolynomial{C}) where C
    return LP1 + (LP2 * (-1))
end

@inline function Base.:-(LP::LabelledPolynomial{C}) where C
    return LabelledPolynomial(LP.index, LP.signature, (LP.poly * (-1)))
end

@inline function Base.:(//)(LP::LabelledPolynomial{C}, c::Number) where C
    return LabelledPolynomial(LP.index, LP.signature, LP.poly * (1//c))
end

@inline Base.iszero(LP::LabelledPolynomial) = iszero(LP.poly)

# --- Comparison ---

# Comparison for sorting lists of LabelledPolynomials
@inline function Base.isless(LP1::LabelledPolynomial, LP2::LabelledPolynomial)
    if LP1.index != LP2.index
        return LP1.index < LP2.index
    end
    return LP1.signature < LP2.signature
end

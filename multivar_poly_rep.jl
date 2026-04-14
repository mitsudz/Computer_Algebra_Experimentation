#=
# Fast representation of polynomials
#
# Note: Gemini 3 Flash used to help generate code
=#

# TODO - documentation (production grade), 
# TODO - error checking for more than 7 variables, debug mode add-ons if needed
# TODO - use @inbounds for array accesses
# TODO - scalar division would be useful for things other than rational (so we can to finite fields later)
# TODO - tests need to be made more comprehensive
# TODO - can num_vars become a value type?


using DynamicPolynomials


# --- Bit Manipulation Helpers ---
const BITS_PER_VAR = 16
const TOTAL_BITS = 128
const MAX_EXP = (1 << BITS_PER_VAR) - 1
const OVERFLOW_MASK = 0x80808080808080808080808080808080 # Every 16th bit set 

# --- Packed Bit Representation for Monomial in Graded Lexical Order (see Maple TODO - REF)

# We'll use UInt128 to store: [Total Degree (16b)][x1 (16b)][x2 (16b)]...[x5 (16b)]
# The last bit of each xi (and total degree) is an overflow bit. Effectively we have a maximum exponent of 2^15 ~ 32K
# Then for testing divisibility we simply need to subtract and check for overflow 
# (we will assume overflow cannot occur on addition)
# This fits 7 variables + Total Degree.
struct GrLexMonomial
    bits::UInt128
end

@inline identity_monomial() = GrLexMonomial(UInt128(0))

@inline total_degree(m::GrLexMonomial) = return Int(m.bits >> (TOTAL_BITS - BITS_PER_VAR))

@inline Base.hash(m::GrLexMonomial, h::UInt) = hash(m.bits, h)


struct Term{C}
    coeff::C
    mono::GrLexMonomial
end

@inline total_degree(t::Term) = total_degree(t.mono)

struct FastPoly{C}
    # Keep terms sorted by mono.bits DESCENDING for GrLex
    terms::Vector{Term{C}}
end

@inline function total_degree(p::FastPoly)
    iszero(p) && return 0
    return total_degree(leading_monomial(p))
end

# --- Getters (Inlined) ---

@inline leading_term(p::FastPoly) = p.terms[1]
@inline leading_monomial(p::FastPoly) = p.terms[1].mono
@inline leading_coefficient(p::FastPoly) = p.terms[1].coeff

# --- Zero Functionality ---

@inline Base.zero(::Type{FastPoly{C}}) where C = FastPoly(Term{C}[])
@inline Base.zero(p::FastPoly{C}) where C = zero(FastPoly{C})
FastPoly{C}() where C = zero(FastPoly{C})


# --- LCM and SIMD-like operations ---

"""
Computes LCM for GrLex packed monomials.
Calculates component-wise max and updates the total degree.
"""
@inline function lcm(m1::GrLexMonomial, m2::GrLexMonomial, num_vars::Int)
    b1 = m1.bits
    b2 = m2.bits

    new_bits = UInt128(0)
    new_total_deg = 0

    # We iterate through the packed slots.
    # For GrLex, we want max(exp_i, exp_j)
    for i in 1:num_vars
        shift = 128 - (i + 1) * BITS_PER_VAR

        # Vectorized-style extraction using masking
        v1 = (b1 >> shift) & MAX_EXP
        v2 = (b2 >> shift) & MAX_EXP

        mv = max(v1, v2)
        new_total_deg += Int(mv)
        new_bits |= (UInt128(mv) << shift)
    end

    # Set the new total degree in the highest 16 bits
    new_bits |= (UInt128(new_total_deg) << (128 - BITS_PER_VAR))

    return GrLexMonomial(new_bits)
end

"""
Subtracts the exponents of m2 from m1.
Equivalent to m1 / m2 in monomial arithmetic.

Precondition: m2 | m1
"""
@inline div_multiple(m1::GrLexMonomial, m2::GrLexMonomial) = GrLexMonomial(m1.bits - m2.bits)

"""
Finds monomials u and v such that u*LM(f) = v*LM(g) = lcm(LM(f), LM(g)).
"""
@inline function critical_pair(f::FastPoly{C}, g::FastPoly{C}, num_vars::Int) where C
    lm_f = leading_monomial(f)
    lm_g = leading_monomial(g)

    lcm_fg = lcm(lm_f, lm_g, num_vars)

    u = div_multiple(lcm_fg, lm_f)
    v = div_multiple(lcm_fg, lm_g)

    return (u, v)
end

"""
Checks whether m1 | m2. 

Assumes m1 and m2 are valid and no overflow has occurred iin the representation already
"""
@inline function divides(m1::GrLexMonomial, m2::GrLexMonomial)
    # TODO - WRITE A PROPER BIT OF DOCUMENTATION HERE EXPLAINING THE BIT HACKING:
    # 1. Subtraction: If m1[i] > m2[i], the borrow bit at the top of the slot flips
    # 2. XOR/OR: Check if any safety bits are active in the result
    # (Actually, a simpler way is to check if (m2 - m1) wiped out any safety bits or set them)

    # Standard 'SIMD-within-a-register' subtraction check:
    # We use a mask of all "low bits" (the safety bits).
    diff = (m2.bits | OVERFLOW_MASK) - m1.bits
    return (diff & OVERFLOW_MASK) == OVERFLOW_MASK
end

# --- Comparison Operations ---

@inline Base.:(==)(m1::GrLexMonomial, m2::GrLexMonomial) = m1.bits == m2.bits
@inline Base.:(==)(t1::Term{C}, t2::Term{C}) where C = t1.coeff == t2.coeff && t1.mono == t2.mono
@inline Base.:(==)(p1::FastPoly{C}, p2::FastPoly{C}) where C = p1.terms == p2.terms
@inline Base.:>(m1::GrLexMonomial, m2::GrLexMonomial) = m1.bits > m2.bits
@inline Base.:<(m1::GrLexMonomial, m2::GrLexMonomial) = m1.bits < m2.bits
@inline Base.:>=(m1::GrLexMonomial, m2::GrLexMonomial) = m1.bits >= m2.bits
@inline Base.:<=(m1::GrLexMonomial, m2::GrLexMonomial) = m1.bits <= m2.bits
@inline Base.isless(m1::GrLexMonomial, m2::GrLexMonomial) = m1 < m2

# --- Core Arithmetic ---

# Addition/Subtraction via a Merge-sort style update
function Base.:+(p1::FastPoly{C}, p2::FastPoly{C}) where C
    new_terms = Term{C}[]
    i, j = 1, 1
    while i <= length(p1.terms) && j <= length(p2.terms)
        m1, m2 = p1.terms[i].mono.bits, p2.terms[j].mono.bits
        if m1 == m2
            c = p1.terms[i].coeff + p2.terms[j].coeff
            if !iszero(c)
                push!(new_terms, Term(c, p1.terms[i].mono))
            end
            i += 1; j += 1
        elseif m1 > m2  # GrLex: higher bits come first
            push!(new_terms, p1.terms[i])
            i += 1
        else
            push!(new_terms, p2.terms[j])
            j += 1
        end
    end
    append!(new_terms, p1.terms[i:end])
    append!(new_terms, p2.terms[j:end])
    return FastPoly(new_terms)
end

@inline Base.:-(p1::FastPoly{C}, p2::FastPoly{C}) where C = p1 + (-p2)
@inline Base.:-(p::FastPoly{C}) where C = (-1) * p

@inline function Base.:+(p::FastPoly{C}, t::Term{C}) where C
    iszero(t) && return p
    return p + to_poly(t)
end
@inline Base.:+(t::Term{C}, p::FastPoly{C}) where C = p + t
@inline Base.:-(t::Term{C}) where C = Term{C}(-t.coeff, t.mono)

@inline function Base.:-(p::FastPoly{C}, t::Term{C}) where C
    return p + (-t) 
end

@inline Base.:*(m1::GrLexMonomial, m2::GrLexMonomial) = GrLexMonomial(m1.bits + m2.bits)

# Vectorized multiplication: scalar * monomial * poly
@inline function Base.:*(p::FastPoly{C1}, t::Term{C2}) where {C1, C2}
    PROMOTED_C = promote_type(C1, C2)
    if iszero(t.coeff) 
        return FastPoly(Term{PROMOTED_C}[]) 
    end

    new_terms = Term{PROMOTED_C}[
        Term{PROMOTED_C}(t.coeff * pt.coeff, t.mono * pt.mono)
        for pt in p.terms
    ]
    
    return FastPoly(new_terms)
end
@inline Base.:*(c::Term{C1}, p::FastPoly{C2}) where {C1, C2} = p * c

@inline function Base.:*(p::FastPoly{C}, m::GrLexMonomial) where C
    new_terms = Term{C}[
        Term{C}(pt.coeff, m * pt.mono)
        for pt in p.terms
    ]
    return FastPoly(new_terms) 
end
@inline Base.:*(c::GrLexMonomial, p::FastPoly{C}) where C = p * c

# Vectorized multiplication: FastPoly * Scalar
@inline function Base.:*(p::FastPoly{C1}, c::C2) where {C1, C2}
    PROMOTED_C = promote_type(C1, C2)
    iszero(c) && return FastPoly(Term{PROMOTED_C}[])

    new_terms = map(p.terms) do pt
        Term{PROMOTED_C}(pt.coeff * c, pt.mono)
    end

    return FastPoly{PROMOTED_C}(new_terms)
end
@inline Base.:*(c::Number, p::FastPoly{C}) where C = p * c

@inline Base.://(p::FastPoly{C1}, c::C2) where {C1, C2} = p * (1//c)

@inline Base.iszero(p::FastPoly{C}) where C = isempty(p.terms)
@inline Base.iszero(t::Term{C}) where C = iszero(t.coeff)

# --- Conversion ---

function Base.convert(::Type{FastPoly{C}}, m::GrLexMonomial) where C
    return FastPoly([Term(one(C), m)])
end

@inline to_poly(t::Term{C}) where C = FastPoly([t])

@inline to_poly(m::GrLexMonomial, coeff::C = 1) where C = FastPoly([Term(coeff, m)])

function pack_exponents(exps::Vector{Int})
    total_deg = sum(exps)
    # Start with total degree in the highest 16 bits
    # Shifted to the very top of the 128-bit word
    packed = UInt128(total_deg) << (128 - BITS_PER_VAR)
    
    for (i, e) in enumerate(exps)
        shift = 128 - (i + 1) * BITS_PER_VAR
        packed |= (UInt128(e) << shift)
    end
    return GrLexMonomial(packed)
end

function unpack_exponents(m::GrLexMonomial, num_vars::Int)
    exps = zeros(Int, num_vars)
    for i in 1:num_vars
        shift = 128 - (i + 1) * BITS_PER_VAR
        exps[i] = Int((m.bits >> shift) & MAX_EXP)
    end
    return exps
end

function from_dynamic(poly, vars)
    v_map = Dict(v => i for (i, v) in enumerate(vars))
    n = length(vars)
    @assert n <= 7 "Error: FastPoly does not yet support more than 7 indeterminants"

    # 1. Detect the coefficient type (e.g., Int64 or Rational{BigInt})
    C = coefficient_type(poly)

    # 2. Initialize the vector with the CONCRETE type C
    t_list = Term{C}[]

    for t in terms(poly)
        exps = zeros(Int, n)
        for (v, e) in powers(monomial(t))
            exps[v_map[v]] = e
        end
        # Now pushing Term{C} into Vector{Term{C}} works perfectly
        push!(t_list, Term(C(coefficient(t)), pack_exponents(exps)))
    end

    # Ensure GrLex sort
    sort!(t_list, by=x->x.mono.bits, rev=true)
    return FastPoly(t_list)
end

function to_dynamic(fp::FastPoly{C}, vars) where C
    n = length(vars)
    poly = zero(vars[1]) * zero(C)
    for t in fp.terms
        exps = unpack_exponents(t.mono, n)
        m = prod(vars[i]^exps[i] for i in 1:n)
        poly += t.coeff * m
    end
    return poly
end

function generate_converter(vars) 
    function convert_poly(input)
        _to_fast(p) = from_dynamic(p, vars)
        _to_dynamic(p) = to_dynamic(p, vars)

        if input isa AbstractVector && !(input isa DynamicPolynomials.AbstractPolynomial)
            if isempty(input) 
                return input 
            end
            return first(input) isa FastPoly ? [_to_dynamic(p) for p in input] : [_to_fast(p) for p in input]

        elseif input isa FastPoly
            return _to_dynamic(input)

        elseif input isa DynamicPolynomials.AbstractPolynomial
            return _to_fast(input)

        else
            error("Unsupported input type.")
        end
    end

    return convert_poly
end


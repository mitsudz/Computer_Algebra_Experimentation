#
# First implementation of F5B algorithm (Sun, Wang, 2010/11)
#

# TODO - REPLACE ALL THIS CRAP WITH MY OWN (FAST) UNDERLYING POLYNOMIAL REPRESENTATION

include("grobner_helpers.jl")
import MultivariatePolynomials: leading_monomial, leading_term, divides
using DataStructures
import Base:+, -, //

struct LabelledPolynomial{P<:AbstractPolynomialLike, T<:AbstractMonomial}
    index::Int
    signature::T
    poly::P
end

function is_sig_greater(LP1::LabelledPolynomial, LP2::LabelledPolynomial)
    if LP1.index == LP2.index
        return LP1.signature > LP2.signature
    end
    return LP1.index < LP2.index
end

function LabelledPolynomial(i::Int, sign::T, p::P) where {P, T}
    return LabelledPolynomial{P, T}(i, sign, p)
end

""" Leading monomial of poly(F), F a labelled polynomial """
function leading_monomial(F::LabelledPolynomial{P}) where {P <: AbstractPolynomial}
    return leading_monomial(F.poly)
end # leading_monomial

""" Leading term of poly(F), F a labelled polynomial """
function leading_term(F::LabelledPolynomial{P}) where {P <: AbstractPolynomial}
    return leading_term(F.poly)
end # leading_monomial

#=
"""
Multiplication of a LabelledPolynomial with a monomial.
"""
function mult(LP::LabelledPolynomial{P, T}, p::AbstractPolynomial) where {P <: AbstractPolynomial, T<:AbstractMonomial}
    @assert nterms(p) == 1 "F5B Error: LabelledPolynomial must be multplied by a term or monomial only"

    return LabelledPolynomial(
        LP.index, 
        leading_monomial(p) * LP.signature, 
        leading_term(p) * LP.poly       
    )
end # (*)

function mult(p::AbstractPolynomial, LP::LabelledPolynomial{P, T}) where {P <: AbstractPolynomial, T<:AbstractMonomial}
    return mult(LP, p)
end # (*)

"""
Multiplication of a LabelledPolynomial with a constant
"""
function mult(c::Number, LP::LabelledPolynomial)
    return LabelledPolynomial(LP.index, LP.signature, c * LP.poly)
end

function mult(LP::LabelledPolynomial, c::Number)
    return  c * LP
end=#

"""
Multiplication of a LabelledPolynomial with something else. Disregarding types temporarily.
"""
function mult(a::LabelledPolynomial, b)
    return _mult_logic(a, b)
end
function mult(a, b::LabelledPolynomial)
    return _mult_logic(b, a)
end

function _mult_logic(LP::LabelledPolynomial, p)
    if p isa Number
        return LabelledPolynomial(LP.index, LP.signature, p * LP.poly)
    else
        # p is a Monomial, Term, or Polynomial
        return LabelledPolynomial(
            LP.index,
            leading_monomial(p) * LP.signature,
            p * LP.poly
        )
    end
end

"""
Addition of two LabelledPolynomial's.
"""
function (+)(LP1::LabelledPolynomial{P}, LP2::LabelledPolynomial{P}) where {P <: AbstractPolynomial}
    if !is_sig_greater(LP1, LP2)
        return LP2 + LP1
    end
    return LabelledPolynomial(LP1.index, LP1.signature, LP1.poly + LP2.poly)
end # (+)

function (//)(LP::LabelledPolynomial, c)
    !(c isa Number) && @assert false "TRYING TO DIVIDE BY NOT A NUMBER"
    return LabelledPolynomial(LP.index, LP.signature, LP.poly * (1//c))
end

"""
Subtraction of two LabelledPolynomial's.
"""
function (-)(LP1::LabelledPolynomial{P}, LP2::LabelledPolynomial{P}) where {P <: AbstractPolynomial}
    return LP1 + (-LP2)
end # (-)

"""
Comparision of two labelled polynomials.
"""
function isless(LP1::LabelledPolynomial{P}, LP2::LabelledPolynomial{P}) where {P <: AbstractPolynomial}
    throw("Not implemented yet")
end

"""
Negation of a LabelledPolynomial's.
"""
function (-)(LP::LabelledPolynomial{P}) where {P <: AbstractPolynomial}
    return LabelledPolynomial(LP.index, LP.signature, -LP.poly)
end # (-)

# TODO - Testing for correctness, Replace underlying representation, implement tail reduction, switch to mod p
"""
Computes a (reduced) Grobner basis utilising the F5B algorithm.
NOTE: This does not utilise interreduction (see F5C paper for this).
NOTE: See (Sun, Wang, 2010/11) for implementation details.
"""
function f5b(initial_basis::Vector{P})::Vector{P} where {P <: AbstractPolynomial}
    initial_basis = deepcopy(initial_basis)
    initial_basis = sort(initial_basis, by=p -> degree(leading_monomial(p)))
    m = length(initial_basis)
    B = [LabelledPolynomial(i, leading_monomial(one(f)), f) for (i, f) in enumerate(initial_basis)]
    B_ = deepcopy(B)
    CP = PriorityQueue{Tuple{Int, Int}, Int}()

    m = length(B)
    for i in 1:m
        for j in i+1:m
            u, v = crit_pair(B[i].poly, B[j].poly)
            # Sugar degree = degree the S-poly would have before reduction
            sugar = max(maxdegree(u) + maxdegree(B[i].poly), maxdegree(v) + maxdegree(B[j].poly))
            enqueue!(CP, (i, j), sugar)
        end
    end

    #while length(CP) != 0
    #println("Initial CP: ", length(CP))
    curr_deg = 0
    printed = false
    while !isempty(CP)
        ((F_gen, G_gen), sugar) = dequeue_pair!(CP)
        F, G = B[F_gen], B[G_gen]

        if curr_deg != sugar
            #println("Processing degree: ", sugar)
            curr_deg = sugar
        end
    
        u, v = crit_pair(F.poly, G.poly)
        uF, vG = mult(u, F), mult(v, G)

        #### SIGNATURE DROP - TODO - INVESTIGATE SIGNATURE DROP FURTHER!!! Def 2.2 Sun/Wang says assume uF > vG
        #=if uF.index == vG.index && uF.signature == vG.signature # SHOULDN'T HAPPEN FOR REGULAR SEQUENCES!
            println("Signature drop at: S(", F_gen, ", ", G_gen, ")")
            println(length(CP))
        end=#
        (uF.index == vG.index && uF.signature == vG.signature) && continue # Signature drop not allowed
        #### SIGNATURE DROP


        #println(length(B), " ", length(CP), " ", syzygy_criterion(uF, vG, B_), " ", rewritten_criterion(uF, vG, B, F_gen, G_gen))
        if !syzygy_criterion(uF, vG, B) && !rewritten_criterion(uF, vG, B, F_gen, G_gen)
            SP = (uF // coefficient(leading_term(F.poly))) - (vG // coefficient(leading_term(G.poly))) # Labelled S-polynomial
            newP = f5b_reduction(SP, length(B) + 1, B)
            push!(B, newP) # Irrespective of whether it is 0
            printed = false
            if iszero(newP.poly) 
                #println("WARNING: Leaky Criteria")
            end
            if !iszero(newP.poly)
                #append!(CP, [(i, length(B)) for i in 1:(length(B) - 1)]) # Add critical pairs
                new_idx = length(B) # B is now 1 element larger
                for i in 1:(new_idx - 1)
                    B[i].poly == 0 && continue
                    u, v = crit_pair(B[i].poly, B[new_idx].poly)
                    sugar = max(maxdegree(u) + maxdegree(B[i].poly), maxdegree(v) + maxdegree(B[new_idx].poly))
                    enqueue!(CP, (i, new_idx), sugar)
                end
                #println("After adding non-zero, length(CP) = ", length(CP))
            end # if
        end # if

        #=if length(B) == 400
            for p in B
                println(p.index, " ", p.signature)
            end
            error("END")
        end=#
        #=if length(B) == 100*m
            G = reduce_gb([Q.poly for Q in B])
            println(length(B), " ", m, " ", length(G), "\n")
            return f5b(G)
        end=#
        if length(B) % 100 == 0 && !printed
            #println("Basis size is $(length(B))")
            printed = true
        end
    end # while

    #println("END: ", length(B))
    return reduce_gb([Q.poly for Q in B])
end # f5b

"""
Find monomials u/v generating the critical pair.
"""
function crit_pair(f::T, g::T) where {T <: AbstractPolynomial}
    lcm_fg = lcm(f, g)

    # TODO - BEWARE DIV MULTIPLE IS DISGUSTING
    #lt_f = leading_term(f)
    #lt_g = leading_term(g)
    
    lm_f = leading_monomial(f)
    lm_g = leading_monomial(g)


    return (div_multiple(lcm_fg, lm_f), div_multiple(lcm_fg, lm_g))
end # crit_pair

"""
Returns whether uF or vG is divisible by B.
"""
function syzygy_criterion(
        uF::T, 
        vG::T, 
        B::Vector{T})::Bool where {
                                   P <: AbstractPolynomial, 
                                   T <: LabelledPolynomial{P}}
    for LP in B
        iszero(LP.poly) && continue
        if divisible(uF, LP) || divisible(vG, LP)
            return true
        end
    end

    return false
end # syzygy_criterion

"""
Returns whether F is divisible by G
"""
function divisible(F::T, G::T)::Bool where {P <: AbstractPolynomial, T <: LabelledPolynomial{P}}
    #return F.index < G.index && divides(leading_monomial(G), leading_monomial(F))
    return F.index < G.index && divides(leading_monomial(G), F.signature)
end # divisible

"""
Determines whether syzygies from this pair will be handled later.
Returns whether uF or vG is rewritable by B.
"""
function rewritten_criterion(
        uF::T, 
        vG::T, 
        B::Vector{T},
        F_gen::Int,
        G_gen::Int)::Bool where {
                                   P <: AbstractPolynomial, 
                                   T <: LabelledPolynomial{P}}
    for (i, LP) in enumerate(B)
        if (rewritable(uF, F_gen, LP, i) || rewritable(vG, G_gen, LP, i))
            return true
        end
    end

    return false
end # rewritten_criterion

"""
Returns whether F is 'rewritable' by G
"""
function rewritable(
        F::T,
        F_gen::Int, 
        G::T, 
        G_gen::Int)::Bool where {
                                 P <: AbstractPolynomial, 
                                 T <: LabelledPolynomial{P}}
    return F_gen < G_gen && G.index == F.index && divides(G.signature, F.signature)
end # rewritable

"""
Performs F5B reduction on a polynomial by B until it cannot be reduced further.

Invariant: F.signature == f5b_reduction(F).signature
"""
function f5b_reduction(F::T, F_gen::Int, B::Vector{T}) where {P <: AbstractPolynomial, T <: LabelledPolynomial{P}}
    for (G_gen, G) in enumerate(B)
        iszero(G.poly) && continue
        !divides(leading_monomial(G), leading_monomial(F)) && continue

        v = div_multiple(leading_monomial(F), leading_monomial(G))
        vG = mult(v, G)

        !is_sig_greater(F, vG) && continue

        crit_failure = false
        for (i, G_prime) in enumerate(B)
            if divisible(vG, G_prime) || rewritable(vG, G_gen, G_prime, i)
                crit_failure = true
                break
            end
        end
        crit_failure && continue

        c = coefficient(leading_term(F.poly)) / coefficient(leading_term(G.poly))

        # TODO Add assertion for safety that signature doesn't change
        newF = F - mult(c, vG)
        @assert newF.signature == F.signature "F5B Reduction Error - Signature Changed"
        return f5b_reduction(newF, F_gen, B)
    end # for

    return F
end # f5b_reduction

#
# Second implementation of F5B algorithm (Sun, Wang, 2010/11).
# Utilises bespoke multivariate polynomial representation
#
# Author: Mittun Sudhahar
# Institution: The University of Queensland
#

# AI Use: Gemini 3 Flash - Code Writing & Debugging Assistance

include("reduce_gb.jl")
include("multivar_poly_rep.jl")
include("labelled_poly.jl")
using DataStructures

# Online TODO Items:
# 1. Add in lazy evaluations (decide which of the next two to try after doing that)
# 2. Add in tail reductions
# 3. Replace PQ and handle sugar & critical pairs

# TODO - Consider using BitIntegers.jl to get more space for more polynomials
# TODO - Make this a module soon to avoid namespace conflicts
# TODO - Benchmark against old f5b and groebner and buchberger
# TODO - Testing for correctness (unit test each part), Replace underlying representation, 
# TODO - Replace PriorityQueue with something faster and more specialised
# TODO - implement tail reduction, switch to mod p
# TODO - Implement normal selection strategy properly (see F4/F5 or other papers - use REFERENCES!!!)
# TODO - Add in speed ups to how much of the basis is passed in for things like the syzygy criterion etc.
# TODO - Figure out a way to do more stuff in-place
# TODO - Ensure vectorisation is happening if possible in all places (e.g. multiplication by a monomial)
# TODO - Make it so polynomials are ensured to be monic before and after function calls as an invariant
# TODO - Decide on a naming convention for functions vs variables
# TODO - Add a debug mode that can be handled via command line variables at compile time

# --- F5B Core Logic --- #
#= # TODO This should be handled OUTSIDE of this file by a wrapper/main file
function f5b(initial_basis::P, vars) where {P <: AbstractPolynomial}
    initial_basis = sort!(reduce_gb(deepcopy(initial_basis)), by = p -> degree(leading_monomial(p)))

    converter = generate_converter(vars)
    fast_basis = converter(initial_basis)
    groebner_basis = f5b(fast_in, length(vars))
    
    return reduce_gb(converter(groebner_basis))
end
=#

"""
Computes a reduced Grobner basis utilising the F5B algorithm and a bespoke polynomial representation.

Modifies given basis in place.

NOTE: This does not utilise interreduction (see F5C paper (Eder/Perry) for this).
NOTE: See (Sun, Wang, 2010/11) for implementation details.
"""
function _f5b(fast_initial::Vector{FastPoly{C}}, num_vars::Int)::Vector{FastPoly{C}} where C
    #fast_inital = _reduce_gb(fast_initial, num_vars) # occurs in-place # TODO - the new reduce requires a grobner basis
    m = length(fast_initial)

    # [(ei, fi) | i = 1,...,m]
    B = [LabelledPolynomial(i, identity_monomial(), f) for (i, f) in enumerate(fast_initial)]

    # Setup Priority Queue by Sugar Degree (Normal Selection Strategy)
    CP = PriorityQueue{Tuple{Int, Int}, Int}()
    for i in 1:m
        for j in i+1:m
            u, v = critical_pair(B[i].poly, B[j].poly, num_vars)
            enqueue!(CP, (i, j), sugar(u, B[i].poly, v, B[j].poly))
        end #for
    end #for

    # Main Loop
    while !isempty(CP)
        (F_idx, G_idx) = dequeue!(CP)
        F, G = B[F_idx], B[G_idx]

        u, v = critical_pair(F.poly, G.poly, num_vars) # TODO - MAKE IT SO THESE AREN'T RECOMPUTED
        uF = F * u # TODO - LAZY - Delay computation of this multiplication until it is actually necessary 
        vG = G * v

        # Avoid Signature Drop - see Sun/Wang 2010/2011 for why this can be skipped
        uF.index == vG.index && uF.signature == vG.signature && continue

        if !syzygy_criterion(uF, vG, B, num_vars) && !rewritten_criterion(uF, vG, B, F_idx, G_idx, num_vars)
            SP = (uF // leading_coefficient(F.poly)) - (vG // leading_coefficient(G.poly)) # S-polynomial

            newP = f5b_reduction(SP, B, num_vars)
            push!(B, newP) # Irrespective of whether the new labelled polynomial is zero

            # Add new pairs
            if !iszero(newP.poly)
                new_idx = length(B)
                for i in 1:(new_idx - 1)
                    iszero(B[i].poly) && continue

                    u_m, v_m = critical_pair(B[i].poly, B[new_idx].poly, num_vars)
                    enqueue!(CP, (i, new_idx), sugar(u_m, B[i].poly, v_m, newP.poly))
                end #for
            end #if
        end #if
    end #while

    # Return the set of polynomials that form the Basis
    return _reduce_gb([Q.poly for Q in B], num_vars)
end # _f5b

"""
Returns whether uF or vG is divisible by B.
"""
function syzygy_criterion(
        uF::LabelledPolynomial{C},
        vG::LabelledPolynomial{C},
        B::Vector{LabelledPolynomial{C}},
        num_vars::Int)::Bool where {C}
    for LP in B
        iszero(LP.poly) && continue
        if divisible(uF, LP, num_vars) || divisible(vG, LP, num_vars)
            return true
        end #if
    end #for

    return false
end # syzygy_criterion

"""
Returns whether F is divisible by G
"""
@inline function divisible(F::LabelledPolynomial{C}, G::LabelledPolynomial{C}, num_vars::Int)::Bool where C
    return F.index < G.index && divides(leading_monomial(G), F.signature, num_vars)
end # divisible

"""
Determines whether syzygies from this pair will be handled later.
Returns whether uF or vG is rewritable by B.
"""
function rewritten_criterion(
        uF::LabelledPolynomial{C}, 
        vG::LabelledPolynomial{C}, 
        B::Vector{LabelledPolynomial{C}},
        F_gen::Int,
        G_gen::Int,
        num_vars::Int)::Bool where C
    for (i, LP) in enumerate(B)
        if (rewritable(uF, F_gen, LP, i, num_vars) || rewritable(vG, G_gen, LP, i, num_vars))
            return true
        end
    end

    return false
end # rewritten_criterion

"""
Returns whether F is 'rewritable' by G
"""
@inline function rewritable(
        F::LabelledPolynomial{C},
        F_gen::Int, 
        G::LabelledPolynomial{C}, 
        G_gen::Int,
        num_vars::Int)::Bool where C
    return F_gen < G_gen && G.index == F.index && divides(G.signature, F.signature, num_vars)
end # rewritable

"""
Performs F5B reduction on a polynomial by B until it cannot be reduced further.
Assumes we are doing reduction on a potential new polynomial to add to basis B

NOTE: This operation mutates F (or at least I want it to in future)
Invariant: F.signature == f5b_reduction(F).signature
"""
function f5b_reduction( # TODO - Make this whole function mutate in-place!
        F::LabelledPolynomial{C},
        B::Vector{LabelledPolynomial{C}},
        num_vars::Int) where C
    iszero(F.poly) && return F # Early exit
    
    while !iszero(F.poly)
        reduced = false

        for (G_idx, G) in enumerate(B)
            iszero(G.poly) && continue
            !divides(leading_monomial(G), leading_monomial(F), num_vars) && continue

            v = div_multiple(leading_monomial(F), leading_monomial(G)) 
            vG = v * G

            !is_sig_greater(F, vG) && continue

            crit_failure = false
            for (H_idx, H) in enumerate(B)
                if (
                    (!iszero(H.poly) && divisible(vG, H, num_vars)) || 
                    rewritable(vG, G_idx, H, H_idx, num_vars)
                   )
                    crit_failure = true
                    break
                end #if
            end #for
            crit_failure && continue

            # TODO - Make it so all polynomials are monic by the time the enter the basis
            c = leading_coefficient(F.poly) / leading_coefficient(G.poly)

            # TODO - IN-PLACE SUBTRACTION!
            F = F - (c * vG)
            reduced = true
            break
            # @assert newF.signature == F.signature "F5B Reduction Error - Signature Changed" # Only in debug mode
        end #for

        !reduced && break
    end #while

    return F
end # f5b_reduction


# --- Critical Pair Queue Logic --- #

# TODO - Check whether sugar is working correctly - I think it's supposed to take the degree of the original basis p?
"""
Get sugar degree of a critical pair (u, F, v, G)
"""
@inline function sugar(u::GrLexMonomial, F::FastPoly{C}, v::GrLexMonomial,  G::FastPoly{C}) where C
    max(total_degree(u) + total_degree(F), total_degree(v) + total_degree(G))
end


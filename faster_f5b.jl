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

const DEBUG = false

# Online TODO Items:
# 2. Replace PQ and handle sugar & critical pairs
# 3. Add in tail reductions
# 4. Add in a way to do initial inter-reduction (this was a 6x speedup for c5 not counting inter-reduction overhead)

# TODO - Consider using BitIntegers.jl to get more space for more polynomials
# TODO - Consider applying criterion BEFORE adding to the queue - if we do this in the right way, we should only
#        need to test the syzygy criterion when deciding to add a critical pair, and test rewritten when deciding 
#        whether to reduce a critical pair.
# TODO - Make this a module soon to avoid namespace conflicts
# TODO - Note the existence of GVW as a competitor (and winner in certain cases) to F5
# TODO - Benchmark against old f5b and groebner and buchberger
# TODO - Testing for correctness (unit test each part), Replace underlying representation, 
# TODO - Replace PriorityQueue with something faster and more specialised
# TODO - implement tail reduction, switch to mod p
# TODO - Implement normal selection strategy properly (see F4/F5 or other papers - use REFERENCES!!!) 
#        "A New Efficient Algorithm for Computing Gröbner Bases" by Giovini, Mora, Niesi, Robbiano, and Traverso (1991).
# TODO - use sizehint! on SyzygyPool
# TODO - Add in speed ups to how much of the basis is passed in for things like the syzygy criterion etc.
# TODO - Figure out a way to do more stuff in-place
# TODO - Ensure vectorisation is happening if possible in all places (e.g. multiplication by a monomial)
# TODO - Make it so polynomials are ensured to be monic before and after function calls as an invariant
# TODO - Decide on a naming convention for functions vs variables
# TODO - why do I get two leaky criterion on cyclic 5?
# TODO - Add a debug mode that can be handled via command line variables at compile time

# --- Critical Pair Queue Logic --- #

# TODO - Yuck naming
struct CriticalPairQueueElem
    sugar::Int
    i::Int
    j::Int
end
Base.isless(a::CriticalPairQueueElem, b::CriticalPairQueueElem) = a.sugar < b.sugar

# TODO - Check whether sugar is working correctly - I think it's supposed to take the degree of the original basis p?
#        "A New Efficient Algorithm for Computing Gröbner Bases" by Giovini, Mora, Niesi, Robbiano, and Traverso (1991).
"""
Get sugar degree of a critical pair (u, F, v, G)
"""
@inline function sugar(u::GrLexMonomial, F::FastPoly{C}, v::GrLexMonomial,  G::FastPoly{C}) where C
    max(total_degree(u) + total_degree(F), total_degree(v) + total_degree(G))
end


# --- Syzygy Criterion Handling --- #

# TODO - test new syzygy criterion and integrate in and verify correctness,
#        check benchmark speed w/wo initial inter-reduction
#        use the cool new divides operation in way fewer clock cycles
#        try the loop unrolling with Val{N} if possible (otherwise just move to the priority queue etc.)
#        this stuff is too in-depth to start with I think

struct SyzygyPool
    lmonoms::Vector{GrLexMonomial} # This should compile to Vector{UInt128} TODO - @assert isbitstype(sigs) - should keep cache lines hot
    indices::Vector{UInt16} # We assume we never have more than 2^16 initial generators

    function SyzygyPool(l, i)
        @assert isbitstype(GrLexMonomial) "GrLexMonomial must be a bitstype for performance!"
        new(l, i)
    end # SyzygyPool
end

@inline function update_syz_pool(syzygies::SyzygyPool, lm::GrLexMonomial, idx::UInt16)
    push!(syzygies.lmonoms, lm)
    push!(syzygies.indices, idx)
end

"""
Returns whether the syzygy criterion is satisfied by a pair of signatures

Note - in F5B this returns whether uF or vG is divisible by B,
       however, we only require monomial and syzygy index comparisons
       rather than the entire polynomials and hence this is an optimised version.
"""
@inline function syzygy_criterion(
        syzygies::SyzygyPool,
        F_idx::UInt16, F_sig::GrLexMonomial,
        G_idx::UInt16, G_sig::GrLexMonomial
    )
    lmonoms = syzygies.lmonoms # The leading monomial of some polynomial
    indices = syzygies.indices
    @inbounds for i in eachindex(lmonoms)
    	lm = lmonoms[i]
    
        if F_idx < indices[i] && divides(lm, F_sig) 
            return true
    	end #if

        if G_idx < indices[i] && divides(lm, G_sig)
            return true
    	end #if
	end #for
    return false
end # syzygy_criterion

"""
Returns whether the syzygy criterion is satisfied by a single signatures

Note - in F5B this returns whether vG is divisible by B,
       however, we only require monomial and syzygy index comparisons
       rather than the entire polynomials and hence this is an optimised version.
"""
@inline function syzygy_criterion(
        syzygies::SyzygyPool,
        idx::UInt16, sig::GrLexMonomial
    )
    lmonoms = syzygies.lmonoms # The leading monomial of some polynomial
    indices = syzygies.indices
    @inbounds for i in eachindex(lmonoms)
        lm = lmonoms[i]
        if idx < indices[i] && divides(lm, sig) 
            return true
    	end #if
	end #for
    return false
end # syzygy_criterion


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
    #fast_inital = _reduce_gb(fast_initial) # occurs in-place # TODO - the new reduce requires a grobner basis
    m = length(fast_initial)

    # [(ei, fi) | i = 1,...,m]
    B = [LabelledPolynomial(i, identity_monomial(), f) for (i, f) in enumerate(fast_initial)]

    # Setup Priority Queue by Sugar Degree (Normal Selection Strategy)
    CP = BinaryHeap{CriticalPairQueueElem}(Base.Order.Forward) # Min-Heap TODO - replace with bucket queue
    for i in 1:m
        for j in i+1:m
            u, v = critical_pair(B[i].poly, B[j].poly, num_vars)
            push!(CP, CriticalPairQueueElem(sugar(u, B[i].poly, v, B[j].poly), i, j))
        end #for
    end #for

    # Set up for Syzygy Criterion
    lms = GrLexMonomial[leading_monomial(LP.poly) for LP in B]
	idxs = UInt16[UInt16(LP.index) for LP in B]
	syzygies = SyzygyPool(lms, idxs)

    # Main Loop
    while !isempty(CP)
        cp = pop!(CP)
        F_idx, G_idx = cp.i, cp.j
        F, G = B[F_idx], B[G_idx]

        u, v = critical_pair(F.poly, G.poly, num_vars) # TODO - MAKE IT SO THESE AREN'T RECOMPUTED
        uF_sig = F.signature * u
        vG_sig = G.signature * v

        # Avoid Signature Drop - see Sun/Wang 2010/2011 for why this can be skipped
        F.index == G.index && uF_sig == vG_sig && continue


        if ( !syzygy_criterion(syzygies,
                               UInt16(F.index), uF_sig,
                               UInt16(G.index), vG_sig
                              ) && 
            !rewritten_criterion(uF_sig, F.index, vG_sig, G.index, B, F_idx, G_idx) )
            SP = ((F*u) // leading_coefficient(F.poly)) - ((G*v) // leading_coefficient(G.poly)) # S-polynomial

            newP = f5b_reduction(SP, B, syzygies)
            push!(B, newP) # Irrespective of whether the new labelled polynomial is zero
            #(DEBUG && iszero(newP.poly)) && println("WARNING: leaky criterion") 

            # Add new pairs
            if !iszero(newP.poly)
                update_syz_pool(syzygies, leading_monomial(newP.poly), UInt16(newP.index))

                new_idx = length(B)
                for i in 1:(new_idx - 1)
                    iszero(B[i].poly) && continue

                    u_m, v_m = critical_pair(B[i].poly, B[new_idx].poly, num_vars)
                    push!(CP, CriticalPairQueueElem(sugar(u_m, B[i].poly, v_m, newP.poly), i, new_idx))
                end #for
            end #if
        end #if
    end #while

    # Return the set of polynomials that form the Basis
    return _reduce_gb([Q.poly for Q in B])
end # _f5b


"""
Determines whether syzygies from this pair will be handled later.
Returns whether uF or vG is rewritable by B.
"""
function rewritten_criterion(
        #uF::LabelledPolynomial{C}, 
        #vG::LabelledPolynomial{C}, 
        uF_sig::GrLexMonomial,
        F_idx::Int,
        vG_sig::GrLexMonomial,
        G_idx::Int,
        B::Vector{LabelledPolynomial{C}},
        F_gen::Int,
        G_gen::Int)::Bool where C
    for (i, LP) in enumerate(B) # TODO - we don't need to iterate through all of them this way since we only need things generated later!
        LP_sig = LP.signature
        LP_idx = LP.index
        if (rewritable(uF_sig, F_idx, F_gen, LP_sig, LP_idx, i) || 
            rewritable(vG_sig, G_idx, G_gen, LP_sig, LP_idx, i))
            return true
        end
    end

    return false
end # rewritten_criterion

"""
Returns whether F is 'rewritable' by G
"""
@inline function rewritable(
        #F::LabelledPolynomial{C},
        F_sig::GrLexMonomial,
        F_idx::Int,
        F_gen::Int, 
        #G::LabelledPolynomial{C}, 
        G_sig::GrLexMonomial,
        G_idx::Int,
        G_gen::Int)::Bool
    return F_gen < G_gen && G_idx == F_idx && divides(G_sig, F_sig)
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
        syzygies::SyzygyPool) where C
    iszero(F.poly) && return F # Early exit
    
    while !iszero(F.poly)
        reduced = false

        for (G_idx, G) in enumerate(B)
            iszero(G.poly) && continue
            !divides(leading_monomial(G), leading_monomial(F)) && continue

            v = div_multiple(leading_monomial(F), leading_monomial(G)) 
            #vG = v * G # TODO - DELAY COMPUTATION!
            vG_sig = v * G.signature

            !is_sig_greater(F.signature, F.index, vG_sig, G.index) && continue

            syzygy_criterion(syzygies, UInt16(G.index), vG_sig) && continue

            crit_failure = false
            for (H_idx, H) in enumerate(B)
                if (
                    rewritable(vG_sig, G.index, G_idx, H.signature, H.index, H_idx)
                   )
                    crit_failure = true
                    break
                end #if
            end #for
            crit_failure && continue

            # TODO - Make it so all polynomials are monic by the time the enter the basis
            c = leading_coefficient(F.poly) / leading_coefficient(G.poly)

            # TODO - IN-PLACE SUBTRACTION!
            F = F - (c * (v * G))
            reduced = true
            break
            # @assert newF.signature == F.signature "F5B Reduction Error - Signature Changed" # Only in debug mode
        end #for

        !reduced && break
    end #while

    return F
end # f5b_reduction



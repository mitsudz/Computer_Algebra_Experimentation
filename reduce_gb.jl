# 
# File for helper functions related to computing Grobner basis
# Implements functionality based on bespoke multivariate late polynomial representations
#

include("multivar_poly_rep.jl")

# TODO - Make my helper functions start with _ or __ to avoid namespace issues (alongside making this a module)

"""
Computes a remainder after multivariate polynomial division of `f` by the `divisors`.
"""
function _polyrem(f::FastPoly{C}, divisors::Vector{FastPoly{C}}, num_vars::Int) where C
    r = zero(f)
    while !iszero(f)
        lm_f = leading_monomial(f)
        found_divisor = false

        # Divide leading term by a divisor and update f
        for divisor in divisors
            lm_div = leading_monomial(divisor)
            if divides(lm_div, lm_f, num_vars)
                found_divisor = true
                coeff_ratio = leading_coefficient(f) // leading_coefficient(divisor)
                f -= divisor * div_multiple(lm_f, lm_div) * coeff_ratio # TODO - make this in place
                break
            end #if
        end #for

        # Add to remainder if no divisor found
        if found_divisor == false
            r += leading_term(f)
            f -= leading_term(f) # TODO - make these both in-place - THIS IS BAD HERE ESPECIALLY!
        end #if
    end #while

    return r
end # polyrem

"""
Computes reduced gröbner basis. 
Modifies input basis G in place.
Returns reduced Grobner basis in ascending order of leading monomial.

Precondition: G must be a Gröbner basis
"""
function _reduce_gb(G::Vector{FastPoly{C}}, num_vars::Int)::Vector{FastPoly{C}} where C
	G = filter!(p -> !iszero(p), G) # Remove 0 polynomials
    isempty(G) && return G # Short circuit edge case

    # Create a minimal Gröbner basis
    # (CLO) C2 S7 Lemma 3: LT(p) in LT<G/{p}> => G/{p} is still Gröbner.
    sort!(G, by = p -> leading_monomial(p))
	tempG = FastPoly{C}[]
    for g in G
    	lm_g = leading_monomial(g)
    	redundant = false
    
    	for p in tempG 
        	if divides(leading_monomial(p), lm_g, num_vars)
            	redundant = true
            	break
        	end #if
    	end #for
    
    	!redundant && push!(tempG, g)
	end #for
	G = tempG

    # Create a reduced Gröbner basis
    for i in 1:length(G)
        G[i], G[end] = G[end], G[i]
        G[end] = _polyrem(G[end], G[1:end-1], num_vars) # tail reduce
        G[end] = G[end] // leading_coefficient(G[end]) # make monic
        G[i], G[end] = G[end], G[i]
    end #for

	return G
end # reduce_gb


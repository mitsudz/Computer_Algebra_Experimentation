# 
# File for helper functions related to computing Grobner basis
# Implements functionality based on bespoke multivariate late polynomial representations
#

include("multivar_poly_rep.jl")

# TODO - Make my helper functions start with _ or __ to avoid namespace issues (alongside making this a module)

"""
Computes a remainder after multivariate polynomial division of `f` by the `divisors`.
"""
function polyrem(f::FastPoly{C}, divisors::Vector{FastPoly{C}}, num_vars::Int) where C
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
Computes reduced (gröbner) basis. 
Modifies input basis G in place.
"""
function reduce_gb(G::Vector{FastPoly{C}}, num_vars::Int)::Vector{FastPoly{C}} where C
	G = filter!(p -> !iszero(p), G) # Remove 0 polynomials

	# Make all monic
	map!( p -> p // leading_coefficient(p), G, G)
    length(G) == 1 && return G # Short circuit edge case TODO - does length(G) == 0 matter?

    sort!(G, by = p -> leading_monomial(p))

	# Repeatedly reduce G until it is stable
	while true
		changed = false

        # gi -> gi ÷ (G - {gi})
		i = 1
		while i <= length(G)
			G[i], G[end] = G[end], G[i]
			temp = polyrem(G[end], G[1:end-1], num_vars)
			if temp != G[end]
				changed = true
				G[end] = temp
			end #if

			# Remove redundant polynomial
			if iszero(G[end])
				changed = true
				pop!(G)
				continue
			end #if

            # Normalise
			G[end] = G[end] // leading_coefficient(G[end])
			removed = false
			for j in 1:length(G)-1
				if divides(leading_monomial(G[j]), leading_monomial(G[end]), num_vars)
					changed = true
					removed = true
					pop!(G)
					break
				end #if
			end #for
			removed && continue
			
			# Swap back if we are keeping it
			G[i], G[end] = G[end], G[i]
			i += 1
		end #while

		
		!changed && break
	end #while

	sort!(G, by=p -> leading_monomial(p))
	return G
end # reduce_gb


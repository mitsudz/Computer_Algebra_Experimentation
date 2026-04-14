# File for helper functions related to computing Grobner basis

import Base: lcm
using MultivariatePolynomials

# NOTE - For simplicity I am currently assuming we are working over a field -> currently I'm using the rationals

# Note - in future use @inbounds to optimise array accesses, and figure out how to parallelise this algorithm
# https://juliaalgebra.github.io/MultivariatePolynomials.jl is my new source of truth for this package:)

"""
Computes the least common multiple of the leading terms of two polynomials.
"""
function lcm(f::T, g::T)::T where {T <: AbstractPolynomial}
    lm_f = MultivariatePolynomials.leading_monomial(f)
    lm_g = MultivariatePolynomials.leading_monomial(g)
    lcm(lm_f, lm_g)
end

"""
Computes the S-polynomial of `f` and `g`
spoly(f, g) = (lcm/lt(f))*f - (lcm/lt(g))*g
"""
function spoly(f::T, g::T)::T where {T <: AbstractPolynomial}
    lcm_fg = lcm(f, g)

    lt_f = MultivariatePolynomials.leading_term(f)
    lt_g = MultivariatePolynomials.leading_term(g)

    return MultivariatePolynomials.div_multiple(lcm_fg, lt_f) * f - MultivariatePolynomials.div_multiple(lcm_fg, lt_g) * g
end

"""
Computes a remainder after multivariate polynomial division of `f` by the divisors.
"""
function polyrem(f::T, divisors::Vector{T})::T where {T <: AbstractPolynomial}
    r = 0*f
    while f != 0
        lt_f = MultivariatePolynomials.leading_term(f)
        found_divisor = false

        #println("Leading term ", lt_f)
        # Divide leading term by a divisor and update f
        for divisor in divisors
            #println("Divisor ", divisor)
            lt_div = MultivariatePolynomials.leading_term(divisor)
            #println("Divisor leading term: ", lt_div)
            if MultivariatePolynomials.divides(lt_div, lt_f)
				#println("divisor was $(divisor)")
                #println("leading term of divisor $lt_div divides $lt_f")
                found_divisor = true
                #println("f was $f")
                f -= divisor * MultivariatePolynomials.div_multiple(lt_f, lt_div)
                #println("f is now $f")
                break
            end
        end

        # Add to remainder if no divisor found
        if found_divisor == false
            r += lt_f
            f -= lt_f
        end
    end

    return r
end

function reduce_gb(G::Vector{T})::Vector{T} where {T <: AbstractPolynomial}
	# Remove 0 polynomials
	G = filter(p -> !iszero(p), G)

	# Make all monic
	map!( p -> p / MultivariatePolynomials.leading_coefficient(p), G, G)
	length(G) == 1 && return G # Short circuit edge case

	# Repeatedly reduce G until it is stable
	while true
		changed = false

		# gi -> gi ÷ (G / gi)
		i = 1
		while i <= length(G)
			#@show G[i]
			G[i], G[end] = G[end], G[i]
			temp = polyrem(G[end], G[1:end-1])
			if temp != G[end]
				changed = true
				G[end] = temp
			end

			# Remove redundant polynomial
			if iszero(G[end])
				changed = true
				pop!(G)
				continue
			end
			G[end] = G[end] / MultivariatePolynomials.leading_coefficient(G[end])
			removed = false
			for j in 1:length(G)-1
				if MultivariatePolynomials.divides(MultivariatePolynomials.leading_term(G[j]), 
                                                   MultivariatePolynomials.leading_term(G[end]))
					changed = true
					removed = true
					pop!(G)
					break
				end
			end
			removed && continue
			
			# Swap back if we are keeping it
			G[i], G[end] = G[end], G[i]
			i += 1
		end

		
		!changed && break
	end

	sort!(G, by=p -> MultivariatePolynomials.leading_term(p))
	return G
end

function change_coeff_type(::Type{C}, f) where {C}
	sum(C(c)*m for (c,m) in zip(coefficients(f), monomials(f)))
end


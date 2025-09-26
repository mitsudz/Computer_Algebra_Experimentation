#= 
TODO - Write my current understanding of Faugere's F4.
As of right now it seems like it centres for correctness around the theorem
that taking the reduced matrix and combining the row echelon basis with the original
to get all unique leading terms will essentially give us all the information
the S-polynomials we chose could possibly give us.

Apparently these mstrices are likely to be sparse so sparse matrix
reduction algorithms will be useful (makes sense ig given sparsity is more
likely with more variables). Probably just try and call FLINT on this one.
=#

include("grobner_helpers.jl")
include("heap.jl")
using LinearAlgebra
using RowEchelon

"""
My first attempt at an F4-esque Grobner basis construction algorithm.

First selection strategy will be to literally just choose some proportion
of the remaining S-polynomials (with a min of 4 to ensure termination).
Next will be to use the normal selection strategy as mentioned by Faugere.

Based on pseudo-code from Ideals, Varieties & Algorithms by Cox, Little, 
Shea on basic F4 (Chapter 10, Section 3)
"""
function F4_naive(G::Vector{T})::Vector{T} where {T <: AbstractPolynomial}
	G = deepcopy(G)
	B = [(i, j) for i in 1:length(G) for j in i+1:length(G)]

	iter_n = 0
	while !isempty(B) 
		iter_n += 1
		# Horrendous selection strategy
		m = max(4, length(B) ÷ 3) 
		B_ = B[max(1, end-m+1):end]
		resize!(B, length(B) - min(length(B), m))

		# Compute reduced matrix
		L = [G[i] * (div_multiple(lcm(G[i], G[j]), leading_term(G[i]))) for (i, j) in B_]
		append!(L, [G[j] * (div_multiple(lcm(G[i], G[j]), leading_term(G[j]))) for (i, j) in B_]) # Get both halves of the S-polynomials
		M, monoms, lm_M = computeM(L, G)

		#N = lu(M; check=false).U # Turns out that only doing LU decomp is really slow compared to the RREF?
		# There should be sparse matrix versions of these which are much faster, and specialised ones which satisfy the F4 requirements.
		N = rref(M) # Take row reduced echelon form
		row_mask = [any(N[i, :] .!= 0) for i in 1:size(N,1)]
		N = N[row_mask, :]

		# Convert N to polynomials if required by Faugere's theorem
		# There's probably a bit more filtering based on that theorem but idk
		N_ = [sum(monoms.*N[i,:]) for i in 1:size(N)[1]] 
		filter!(p -> all(m -> !divides(m, leading_monomial(p)), lm_M), N_) 
	
		# Add the new polynomials into G
		for p in N_
			n = length(G)
			push!(G, p)
			for i in 1:n
				push!(B, (i, n+1))
			end
		end
	end # END WHILE

	return G
end

"""
Also as described in the textbook in the same Chapter/Section
"""
function computeM(L::Vector{T}, G::Vector{T}) where {T <: AbstractPolynomial}
	H = L # Don't need a copy
	done = Set([leading_monomial(h) for h in H])
	monoms_to_do = Heap!([m for p in H for m in monomials(p) if !(m in done)])
	monoms_in_heap = Set([m for p in H for m in monomials(p) if !(m in done)])

	while !isempty(monoms_to_do)
		mon = pop!(monoms_to_do)
		delete!(monoms_in_heap, mon)
		push!(done, mon)

		for g in G
			if divides(leading_monomial(g), mon) 
				new_p = g * (div_multiple(mon, leading_monomial(g)))
				push!(H, new_p)
				for m in monomials(new_p)
					if !(m in done) && !(m in monoms_in_heap)
						push!(monoms_to_do, m)
						push!(monoms_in_heap, m)
					end
				end
				break
			end
		end
	end

	# Create matrix from H
	monoms = sort!(collect(done), by=t->leading_term(t), rev=true)
	M = hcat([coefficients(p, monoms) for p in H]...)'
	leading_monoms = [leading_monomial(p) for p in H]
	return M, monoms, leading_monoms
end

;

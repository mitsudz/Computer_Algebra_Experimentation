using DynamicPolynomials, Groebner
using GaloisFields
using Combinatorics # For generating combinations with repetitions (monomials)

const F = @GaloisField 101 # Or 2^31 - 1 Mersenne 2147483647 

const BIG = 0
const FINITE = 1

"""
    generate_random_regular(n::Int, degrees::Vector{Int}; type::Int = KATSURA)

Generates a random, dense homogeneous polynomial system of `n` equations in `n` variables 
over the specified coefficient field. Such systems are pseudo-regular/regular with high probability.
"""
function generate_random_regular(n::Int, degrees::Vector{Int}; type::Int = FINITE)
    length(degrees) == n || error("Must provide exactly $(n) degrees for $(n) variables.")

    # 1. Define variables dynamically to match your other functions
    @polyvar x[1:n]
    system = AbstractPolynomial[]
    
    coeff_type = type == BIG ? big(1)//1 : F(1)

    # 2. Loop to build each of the n polynomials
    for k in 1:n
        target_deg = degrees[k]
        poly = zero(x[1])
        
        # Generate all possible combinations of variables that yield total degree 'target_deg'
        # with_replacement_combinations([1,2,...,n], d) gives exactly the homogeneous monomial exponents
        for comb in with_replacement_combinations(1:n, target_deg)
            # Sample a random coefficient from your Galois Field or Q
            if type == FINITE
                # Random element from the Galois Field
                # (Note: rand(F) returns a field element directly)
                c = rand(F)
            else
                # Random rational coefficient for the BIG type pathway
                c = (rand(-50:50) // 1)
            end
            
            # Avoid packing the system with structural zeros
            if !iszero(c)
                poly += c * prod(x[idx] for idx in comb)
            end
        end
        
        # Guard rail: make sure we didn't end up with a dead zero polynomial
        if iszero(poly)
            poly += coeff_type
        end
        
        push!(system, poly)
    end

    coeff_type = type == BIG ? big(1)//1 : F(1)
    return [p * coeff_type for p in system], x
end

# --- How to use it in your benchmark workflow ---

# Suppose we want a 5-variable system matching a degree profile like Katsura-5 or a custom run
# degrees = [d1, d2, d3, d4, d5]
degrees_profile = [2, 2, 2, 3, 3] 

random_sys, vars = generate_random_regular(5, degrees_profile, type=FINITE)

println("Generated a random pseudo-regular system over: ", coefficient_type(random_sys[1]))
println("First polynomial example: ", random_sys[1])

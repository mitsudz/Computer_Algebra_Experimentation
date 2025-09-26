include("grobner_helpers.jl")

# NOTE - For simplicity I am currently assuming we are working over a field -> currently I'm using the rationals

# Note - in future use @inbounds to optimise array accesses, and figure out how to parallelise this algorithm
# https://juliaalgebra.github.io/MultivariatePolynomials.jl is my new source of truth for this package:)

" Computes a (not necessarily reduced) Grobner basis of I = <g | g in `G`>."
function buchberger(G::Vector{T})::Vector{T} where {T <: AbstractPolynomial}
	G = deepcopy(G)
    while true 
        n = length(G)
        
        # Apply Buchberger criterion to extend G from G'
        for i in 1:n-1
            for j in 2:n
                S = spoly(G[i], G[j])
                r = polyrem(S, G)
                if r != 0
                    # Note - slight modification to original algorithm where we divide by the updated G in this iteration
                    push!(G, r)
                end
            end
        end
        
        # Criterion satisfied if no extensions to G
        if length(G) == n
            break
        end
    end

    return G
end


" Computes a list of remainder after multivariate polynomial division of `f` by the divisors."
function polyrem2(f::T, divisors::Vector{T}) where {T <: AbstractPolynomial}
    r = 0*f
    f_ = 0*f
    
    while f != 0
        lt_f = leading_term(f)
        found_divisor = false

        #println("Leading term ", lt_f)
        # Divide leading term by a divisor and update f
        for divisor in divisors
            #println("Divisor ", divisor)
            lt_div = leading_term(divisor)
            #println("Divisor leading term: ", lt_div)
            if divides(lt_div, lt_f)
                #println("leading term of divisor $lt_div divides $lt_f")
                found_divisor = true
                #println("f was $f")
                f -= divisor * div_multiple(lt_f, lt_div)
                #println("f is now $f")
                break
            end
        end

        # Add to remainder if no divisor found
        if found_divisor == false
            r += lt_f
            if f_ == 0
                f_ = f
            end
            f -= lt_f
        end
    end
    
    if r == 0
        return
    else
        return polyrem2(f_, push!(divisors, r))
    end
end

"""
Computes a (not necessarily reduced) Grobner basis of I = <g | g in `G`>. Does not recompute unnecessary 
S-polynomials according to criterion that division will remain 0 under extended basis
"""
function buchberger2(G::Vector{T})::Vector{T} where {T <: AbstractPolynomial}
	G = deepcopy(G)
    iter = 0 # Debugging variable
    startingPoly = 1
    while true && iter < 10
        iter += 1
        n = length(G)
        
        # Apply Buchberger criterion to extend G from G'
        for i in 1:n-1
            for j in startingPoly:n
                S = spoly(G[i], G[j])
                r = polyrem2(S, G)
            end
        end
        
        # Criterion satisfied if no extensions to G
        if length(G) == n
            break
        end

        # New starting point for S-polynomials
        startingPoly = n+1
    end

    return G
end

#=
# Basic Testing:
@polyvar t
@polyvar u
@polyvar x
@polyvar y
@polyvar z

f = (3//1)x^2 + 2x*y
g = (2//1)x - y
h = (1//1)x^3*y^2 + z^2 + 5z^5*x^3
println(f)
println(g)
@polyvar x y z # starting with lex and 3 variables

S = spoly(f, g)
println(S)

# TODO - TEST THE POLYNOMIAL REMAINDER/DIVISION 
println("Last div should be 0: ", polyrem(f, [(1//1)y^2, (2//1)x-y]))

# TODO - TEST BUCHBERGER'S ALGORITHM!
println(buchberger([f, g, h]))
println(buchberger2([f, g, h]))

a = (1//1)x - t - u
b = (1//1)y - t^2 - 2t*u
c = (1//1)z - t^3 - 3t^2*u
gr = buchberger([a, b, c])

println("\n\n")
for pol in gr
    println(pol)
end
=#


;


include("grobner_helpers.jl")
include("faster_f5b.jl")
include("f5b.jl")
include("buchberger.jl")

using DynamicPolynomials, Groebner

@polyvar x y z t u v

f = (3//1)x^2 + 2x*y
g = (2//1)x - y
h = (1//1)x^3*y^2 + z^2 + 5z^5*x^3
f2 = (3//1)t^2 + 2t*u
g2 = (2//1)t - u
h2 = (1//1)t^3*u^2 + v^2 + 5v^5*u^3*x^2 

G = [f, g, h, f2, g2, h2]

println(buchberger(G))
println("")
println(f5b(G))
println("")
println(groebner(G))

println("TEST 1: Result: ", buchberger(G) == f5b(G) && f5b(G) == groebner(G))

using BenchmarkTools
import Random
@btime 1+1

function generate_easy_system(n_vars, degree, seed=42)
    @polyvar x[1:n_vars]
    Random.seed!(seed)

    system = AbstractPolynomial[]
    for i in 1:n_vars
        p = sum(rand(-5:5) * mono for mono in monomials(x, degree))
        # Ensure the system is "easy" by adding a dominant term
        p += 50 * x[i]^degree
        push!(system, p)
    end
    return [p * (big(1)//big(1)) for p in system]
end
#=
for i in 1:10
    G_ = generate_easy_system(2, 3, 42*i)
    println("\n\n\n\n Initial Basis:")
    println(G_)
    println("\n\n\n\n My F5B Result After Reduction:")
    println(f5b(G_))
    println("\n\n\n\n Actual Reduced Groebner Basis:")
    println(groebner(G_))
    println("\n\n\n\n")
    @assert groebner(G_) == f5b(G_) "Test $i Failed"
    println("Test $i Passed")
end
=#


#=
println("")
result = @benchmark buchberger($G)
println("Naive Buchberger")
show(stdout, "text/plain", result)
println("")
result = @benchmark f5b($G)
println("F5B Simplistic Version")
show(stdout, "text/plain", result)
println("")
result = @benchmark groebner($G)
println("Groebner.jl")
show(stdout, "text/plain", result)
println("")=#

"""
Generates the Cyclic-n system.
For 5 & 7 are regular, 4 & 6 are not 
"""
function generate_cyclic(n::Int)
    @polyvar x[1:n]
    system = AbstractPolynomial[]

    # Generate the first n-1 equations
    for k in 1:n-1
        # Sum of products of k consecutive variables (circularly)
        eq = sum(prod(x[mod1(i + j, n)] for j in 0:k-1) for i in 1:n)
        push!(system, eq)
    end

    # The final equation: x1*x2*...*xn - 1
    push!(system, prod(x) - 1)

    return [p * (1//1) for p in system]
end

c5 = generate_cyclic(5)
c5 = [map_coefficients(c -> Rational{BigInt}(c), p) for p in c5]
println("\n\n***** Cyclic 5 *****\n")
println("TEST 2: Result = ", groebner(c5) == f5b(c5))

#=
println("\nNow completing tests on Cyclic 5 (which is regular):\n")
#result = @benchmark buchberger($c5)
#=println("Naive Buchberger")
show(stdout, "text/plain", result)
println("\n")=#
#result = @benchmark f5b($c5)
println("F5B Simplistic Version")
#show(stdout, "text/plain", result)
println("\n")
#result = @benchmark groebner($c5)
println("Groebner.jl")
#show(stdout, "text/plain", result)
println("\n")
println(groebner(c5) == f5b(c5))
=#


function generate_katsura(n)
    # We need n+1 variables: u0 to un
    @polyvar u[1:n+1] # u[1] is u0, u[2] is u1, etc.
    system = AbstractPolynomial[]

    # Linear equation: u0 + 2*u1 + 2*u2... + 2*un - 1 = 0
    linear_eq = u[1] + 2 * sum(u[2:end]) - 1
    push!(system, linear_eq)

    # Remaining n equations
    for m in 0:n-1
        # The formula: sum_{i=-n}^n u_{|i|} * u_{|m-i|} = u_m
        # Because u_i = u_{-i}, we simplify to indices 0...n
        term_sum = zero(u[1])
        for i in -n:n
            idx1 = abs(i) + 1
            idx2 = abs(m - i) + 1
            
            # Only add if the index is within our variable range [0, n]
            if idx1 <= n+1 && idx2 <= n+1
                term_sum += u[idx1] * u[idx2]
            end
        end
        push!(system, term_sum - u[m+1])
    end

    return system
end


# TODO - benchmarkable (to remove setup time), show, $ sign, use the actual tests, 
#        get the actual F5 but use
#        Groebner.jl's matrix reduction technology. Use ProfileView if necessary.
#        Find large regular sequences. Start by upgrading to F5C.
#        Substitute in a better underlying polynomial representation (maybe borrow this 
#        also from Groebner.jl)
#        Re-tag this to the relevant git repo. Get proper docs up and use better tooling

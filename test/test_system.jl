using DynamicPolynomials, Groebner

#=
@polyvar x y z t u v

f = (3//1)x^2 + 2x*y
g = (2//1)x - y
h = (1//1)x^3*y^2 + z^2 + 5z^5*x^3
f2 = (3//1)t^2 + 2t*u
g2 = (2//1)t - u
h2 = (1//1)t^3*u^2 + v^2 + 5v^5*u^3*x^2 

G = [f, g, h, f2, g2, h2]
=#

#=
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

    return [p * (big(1)//1) for p in system], x
end

c5, vars = generate_cyclic(5)
=#

# TODO - Try with Katsura after making some further improvements if we want to go past C5.

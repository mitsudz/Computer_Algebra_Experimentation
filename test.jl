using Groebner, DynamicPolynomials
include("f5b.jl")

@polyvar x y

#f = (1351//1)y^3 - 103x*y^2 + 26x^2*y
#g = (1//1)y^3 - x*y^2 + 52x^3
f = (1//1)y^3 - x*y^2 + x^2*y
g = (1//1)y^3 - x*y^2 + x^3

G = [f, g]

f5b(G)

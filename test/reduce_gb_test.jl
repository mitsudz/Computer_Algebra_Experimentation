using Test

include("../reduce_gb.jl")

@testset "Groebner Helper Tests (FastPoly)" begin
    # Setup variables for a 3-variable system (x, y, z)
    # Using 3 vars to keep bit-masks readable
    num_vars = 3
    
    # Helper to build a monomial quickly
    m(exps) = pack_exponents(vcat(exps, zeros(Int, num_vars - length(exps))))
    
    # Helper to build a poly from list of (coeff, exps)
    p(terms) = FastPoly([Term(Rational{BigInt}(c), m(e)) for (c, e) in terms])

    @testset "Multivariate Division (polyrem)" begin
        # Test 1: Basic Division
        # f = x^2 - y^2, divisors = [x - y] -> rem = 0
        f = p([(1, [2, 0, 0]), (-1, [0, 2, 0])])
        d = [p([(1, [1, 0, 0]), (-1, [0, 1, 0])])]
        @test iszero(polyrem(f, d, num_vars))

        # Test 2: Division with Non-Zero Remainder
        # f = x^2*y + x*y^2 + y^2, divisors = [x^2, y]
        # x^2*y -> 0 (via x^2)
        # x*y^2 -> 0 (via y)
        # y^2   -> 0 (via y)
        f2 = p([(1, [2, 1, 0]), (1, [1, 2, 0]), (1, [0, 2, 0])])
        d2 = [p([(1, [2, 0, 0])]), p([(1, [0, 1, 0])])]
        @test iszero(polyrem(f2, d2, num_vars))

        # Test 3: Remainder that cannot be reduced
        # f = x + y + 1, divisors = [x^2] -> rem = x + y + 1
        f3 = p([(1, [1, 0, 0]), (1, [0, 1, 0]), (1, [0, 0, 0])])
        d3 = [p([(1, [2, 0, 0])])]
        @test polyrem(f3, d3, num_vars) == f3
    end

    @testset "Basis Reduction (reduce_gb)" begin
        # Test 1: Redundant Elements
        # G = {x^2 - y, x, y} -> Reduced GB should be {x, y}
        G1 = [
            p([(1, [2, 0, 0]), (-1, [0, 1, 0])]),
            p([(1, [1, 0, 0])]),
            p([(1, [0, 1, 0])])
        ]
        red1 = reduce_gb(G1, num_vars)
        @test length(red1) == 2
        @test any(poly -> poly == p([(1, [1, 0, 0])]), red1)
        @test any(poly -> poly == p([(1, [0, 1, 0])]), red1)

        # Test 2: Normalization (Monic)
        # G = {2x - 4y} -> Reduced = {x - 2y}
        G2 = [p([(2, [1, 0, 0]), (-4, [0, 1, 0])])]
        red2 = reduce_gb(G2, num_vars)
        @test leading_coefficient(red2[1]) == 1
        @test red2[1].terms[2].coeff == -2

        # Test 3: Standard Example
        # G = {x + y, x - y} -> Reduced = {x, y}
        G3 = [
            p([(1, [1, 0, 0]), (1, [0, 1, 0])]),
            p([(1, [1, 0, 0]), (-1, [0, 1, 0])])
        ]
        red3 = reduce_gb(G3, num_vars)
        @test length(red3) == 2
        @test iszero(polyrem(p([(1, [1, 1, 0])]), red3, num_vars))
    end
end

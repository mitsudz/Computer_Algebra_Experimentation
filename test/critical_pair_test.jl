using Test
using DynamicPolynomials

include("../faster_f5b.jl")

@testset "Sugar Degree Tests" begin
    @polyvar x1 x2 x3
    vars = [x1, x2, x3]
    C = Rational{Int}

    @testset "Standard S-Poly Sugar" begin
        # F = x1^2 + x2 (Degree 2)
        # G = x1 + x2^2 (Degree 2)
        F = from_dynamic(C(1)x1^2 + x2, vars)
        G = from_dynamic(C(1)x1 + x2^2, vars)
        
        # S-poly(F, G) involves:
        # u = LCM(x1^2, x1)/x1^2 = 1 (Degree 0)
        # v = LCM(x1^2, x1)/x1 = x1 (Degree 1)
        u = pack_exponents([0, 0, 0])
        v = pack_exponents([1, 0, 0])
        
        # sugar = max(0 + 2, 1 + 1) = 2
        @test sugar(u, F, v, G) == 2
    end

    @testset "Imbalanced Degrees" begin
        # F = x1^3 (Degree 3)
        # G = x2 (Degree 1)
        F = from_dynamic(C(1)x1^3, vars)
        G = from_dynamic(C(1)x2, vars)
        
        # LCM(x1^3, x2) = x1^3 * x2
        # u = x2 (Deg 1), v = x1^3 (Deg 3)
        u = pack_exponents([0, 1, 0])
        v = pack_exponents([3, 0, 0])
        
        # sugar = max(1 + 3, 3 + 1) = 4
        @test sugar(u, F, v, G) == 4
    end

    @testset "Higher Variable Selection" begin
        # Testing that sugar correctly picks the max when one side is heavier
        # uF path: Deg 2 + Deg 5 = 7
        # vG path: Deg 4 + Deg 2 = 6
        F_poly = from_dynamic(C(1)x1^5, vars)
        G_poly = from_dynamic(C(1)x2^2, vars)
        u = pack_exponents([2, 0, 0])
        v = pack_exponents([0, 4, 0])
        
        @test sugar(u, F_poly, v, G_poly) == 7
    end
end

using Test
using DynamicPolynomials
using Groebner

include("../faster_f5b.jl")

# TODO - Make into a module first so I can add in buchberger's algorithm as a test reference

# TODO - Move to utilities file
function quick_lp(idx, sig_exps, poly_dynamic, vars)
    p = from_dynamic(poly_dynamic, vars)
    sig = pack_exponents(sig_exps)
    return LabelledPolynomial(idx, sig, p)
end

@testset "f5b_reduction Tests" begin
    @polyvar x1 x2
    vars = [x1, x2]
    num_vars = 2
    C = Rational{Int}

    # Setup Basis: g1 = x1 + 1, g2 = x2 + 1
    g1 = quick_lp(1, [0, 0], C(1)x1 + 1, vars)
    g2 = quick_lp(2, [0, 0], C(1)x2 + 1, vars)
    B = [g1, g2]

    @testset "Basic Top Reduction" begin
        # F = x1*x2, signature (1, x2)
        # Should be reduced by g1 (LM=x1) because (1, x2) > x2 * (1, 1) is false? 
        # No, is_sig_greater(F, vG) must be true.
        # Let's make F have a high signature so it can be reduced.
        F = LabelledPolynomial(1, pack_exponents([0, 2]), from_dynamic(C(1)x1*x2, vars))
        
        # v = x2, G = g1. vG signature is (1, x2).
        # F signature is (1, x2^2). (1, x2^2) > (1, x2), so reduction happens.
        reduced_F = f5b_reduction(F, B, num_vars)
        @test !divides(leading_monomial(g1.poly), leading_monomial(reduced_F.poly), num_vars)
    end

    @testset "Signature Constraint (Prevent Drop)" begin
        # g1 has index 1.
        # F has signature (1, 1). 
        # If we try to reduce F by g1, v = 1. v*g1 has signature (1, 1).
        # is_sig_greater(F, vG) will be false (they are equal).
        # Reduction should NOT occur.
        F = LabelledPolynomial(1, pack_exponents([0, 0]), from_dynamic(C(1)x1 + C(1)x2, vars))
        reduced_F = f5b_reduction(F, B, num_vars)
        
        @test reduced_F.poly == F.poly # No change because of signature rule
    end
end


# TODO - TESTS ARE CURRENTLY FAILING DUE TO NAMESPACE CONFLICTS!!!
@testset "Full F5B Algorithm Convergence" begin
    @polyvar x y z
    vars = [x, y, z]
    converter = generate_converter(vars)
    num_vars = 3
    C = Rational{BigInt}

    @testset "System: Katsura-2" begin
        # Traditional test case
        sys = [x + 2y + 2z - 1, x^2 + 2y^2 + 2z^2 - x, 2x*y + 2y*z - y]
        sys = [C(1)p for p in sys]
        fast_sys = [from_dynamic(C(1)p, vars) for p in sys]
        
        # Your F5B
        my_gb = reduce_gb(converter(_f5b(fast_sys, num_vars)))
        
        # Groebner.jl
        official_gb = Groebner.groebner(sys)
        
        @test my_gb == official_gb
    end

    @testset "System: Cyclic-3" begin
        sys = [x + y + z, x*y + y*z + z*x, x*y*z - 1]
        sys = [C(1)p for p in sys]
        fast_sys = [from_dynamic(C(1)p, vars) for p in sys]
        
        my_gb = reduce_gb(converter(_f5b(fast_sys, num_vars)))
        
        # official_gb_1 = buchberger(sys) # TODO
        #check_ideal_equality(my_gb, official_gb_1, converter) # TODO
        
        official_gb = Groebner.groebner(sys)

        @test my_gb == official_gb
    end

    @testset "System: Cyclic-4 (Complexity Bound)" begin
        @polyvar a b c d
        cvars = [a, b, c, d]
        converter = generate_converter(cvars)
        sys = [a+b+c+d, a*b+b*c+c*d+d*a, a*b*c+b*c*d+c*d*a+d*a*b, a*b*c*d-1]
        sys = [C(1)p for p in sys]
        fast_sys = [from_dynamic(C(1)p, cvars) for p in sys]
        
        my_gb = reduce_gb(converter(_f5b(fast_sys, 4)))
        official_gb = Groebner.groebner(sys)
        
        @test my_gb == official_gb
    end
end

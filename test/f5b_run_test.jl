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

    # Initialize Pool for the tests
    lms = [leading_monomial(LP.poly) for LP in B]
    idxs = UInt16[UInt16(LP.index) for LP in B]
    syzygies = SyzygyPool(lms, idxs)

    @testset "Basic Top Reduction" begin
        # F = x1*x2, signature (1, x2^2)
        # v = x2, G = g1. vG signature is (1, x2).
        # F.sig (1, x2^2) > vG.sig (1, x2), so reduction should happen.
        F = LabelledPolynomial(1, pack_exponents([0, 2]), from_dynamic(C(1)x1*x2, vars))
        
        reduced_F = f5b_reduction(F, B, syzygies)
        
        # Verify g1.LM (x1) no longer divides the result's LM
        @test !divides(leading_monomial(g1.poly), leading_monomial(reduced_F.poly))
    end

    @testset "Signature Constraint (Prevent Drop)" begin
        # g1 has index 1.
        # F has signature (1, 1). 
        # If we try to reduce F by g1, v = 1. v*g1 has signature (1, 1).
        # is_sig_greater(F, vG) will be false (they are equal).
        F = LabelledPolynomial(1, pack_exponents([0, 0]), from_dynamic(C(1)x1 + C(1)x2, vars))
        
        reduced_F = f5b_reduction(F, B, syzygies)
        
        @test reduced_F.poly == F.poly # No change because of signature rule
    end

    @testset "Syzygy Criterion Block" begin
        # Let's create a case where a reduction is logically possible by LM, 
        # but the Syzygy Pool says no.
        
        # Setup: Basis g1 has index 1, LM x1.
        # Suppose we have a reducer vG where vG.sig = x1 and vG.index = 0.
        # This sig is divisible by g1.LM (x1) and 0 < 1. 
        # The pool should trigger a 'continue'.
        
        # F = x1*x2, sig = (10, 1) [Very high index to ensure sig_greater passes]
        F = LabelledPolynomial(10, identity_monomial(), from_dynamic(C(1)x1*x2, vars))
        
        # We are reducing by g1 (index 1). v = x2. 
        # vG signature is (1, x2).
        # Is (1, x2) in pool? pool has (x1, index 1). 
        # x1 does NOT divide x2. No syzygy hit here.
        
        # Now let's try reducing by g1 where v = x1.
        # F = x1^2, sig = (10, 1)
        F2 = LabelledPolynomial(10, identity_monomial(), from_dynamic(C(1)x1^2, vars))
        # Reducer vG: v=x1, G=g1. vG signature is (1, x1).
        # Pool check: does pool.LM(x1) divide vG.sig(x1)? Yes.
        # Is vG.index(1) < pool.index(1)? No (1 < 1 is false).
        
        # To force a hit: Let's add a "Principal Syzygy" to the pool manually
        # that specifically blocks index 0 with LM x2.
        test_lms = [leading_monomial(g1.poly), pack_exponents([0, 1])] # g1.LM and x2
        test_idxs = UInt16[1, 5] # g1 and a "future" generator at index 5
        test_pool = SyzygyPool(test_lms, test_idxs)
        
        # F = x1*x2, sig = (10, 1). Reducer vG = x2*g1.
        # vG signature is (1, x2).
        # Pool check: pool element 2 has LM x2 and index 5.
        # vG.sig x2 is divisible by x2. vG.index 1 < 5.
        # HIT! Reduction should be blocked.
        
        reduced_F2 = f5b_reduction(F, B, test_pool)
        @test reduced_F2.poly == F.poly # Should remain unreduced due to syzygy hit
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
        my_gb = converter(_f5b(fast_sys, num_vars))
        
        # Groebner.jl
        official_gb = Groebner.groebner(sys)
        
        @test my_gb == official_gb
    end

    @testset "System: Cyclic-3" begin
        sys = [x + y + z, x*y + y*z + z*x, x*y*z - 1]
        sys = [C(1)p for p in sys]
        fast_sys = [from_dynamic(C(1)p, vars) for p in sys]
        
        my_gb = converter(_f5b(fast_sys, num_vars))
        
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
        
        my_gb = converter(_f5b(fast_sys, 4))
        official_gb = Groebner.groebner(sys)
        
        @test my_gb == official_gb
    end
end

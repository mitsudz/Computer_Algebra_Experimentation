using Test
using DynamicPolynomials

include("../multivar_poly_rep.jl") 
include("../labelled_poly.jl")

@testset "LabelledPolynomial Robustness Tests (Field: Q)" begin
    @polyvar x1 x2 x3 x4 x5
    vars = [x1, x2, x3, x4, x5]
    
    C = Rational{Int64}
    
    # lp1: sig = x1, poly = 2x1^2 + 3x2
    m_sig1 = pack_exponents([1, 0, 0, 0, 0]) 
    fp1 = from_dynamic(C(2)x1^2 + C(3)x2, vars)
    lp1 = LabelledPolynomial(1, m_sig1, fp1)

    # lp2: sig = x2, poly = x1*x2 + 5
    m_sig2 = pack_exponents([0, 1, 0, 0, 0]) 
    fp2 = from_dynamic(C(1)x1*x2 + C(5), vars)
    lp2 = LabelledPolynomial(1, m_sig2, fp2)

    @testset "Signature Comparison (is_sig_greater)" begin
        @test is_sig_greater(lp1, lp2) == true
        lp_idx2 = LabelledPolynomial(2, m_sig1, fp1)
        @test is_sig_greater(lp1, lp_idx2) == true
    end

    @testset "Monomial and Term Multiplication" begin
        m_x2 = pack_exponents([0, 1, 0, 0, 0])
        res_m = lp1 * m_x2
        
        @test to_dynamic(res_m.poly, vars) == C(2)x1^2*x2 + C(3)x2^2
        @test unpack_exponents(res_m.signature, 5) == [1, 1, 0, 0, 0]
    end

    @testset "Addition and Signature Dominance" begin
        sum_lp = lp1 + lp2
        @test sum_lp.signature.bits == lp1.signature.bits
        @test to_dynamic(sum_lp.poly, vars) == (C(2)x1^2 + C(1)x1*x2 + C(3)x2 + C(5))
        
        # --- FIXED CANCELLATION TEST ---
        # We want to subtract lp2 from something with a HIGHER signature (lp1)
        # We'll make a new poly with sig1 that partially cancels lp2
        # fp_partial: x1^2 - x1*x2 (The x1*x2 will cancel, the x1^2 will remain)
        fp_partial = from_dynamic(C(1)x1^2 - C(1)x1*x2, vars)
        lp_high = LabelledPolynomial(1, m_sig1, fp_partial)
        
        res_cancel = lp_high + lp2
        
        # Should NOT be zero because x1^2 and 5 remain
        @test iszero(res_cancel.poly) == false
        @test to_dynamic(res_cancel.poly, vars) == C(1)x1^2 + C(5)
        @test res_cancel.signature.bits == m_sig1.bits
    end

    @testset "Scalar Division (//)" begin
        # lp1 is 2x1^2 + 3x2. 
        # lp1 // 2 should be 1x1^2 + (3//2)x2
        lp_div = lp1 // 2
        @test leading_coefficient(lp_div.poly) == 1 // 1
        @test lp_div.poly.terms[2].coeff == 3 // 2
        @test lp_div.signature.bits == lp1.signature.bits
    end

    @testset "Type Stability (Rational{BigInt})" begin
        CB = Rational{BigInt}
        fp_rat = from_dynamic(CB(1, 3)*x1, vars)
        lp_rat = LabelledPolynomial(1, m_sig1, fp_rat)
        
        res = lp_rat * big(3)
        @test res.poly.terms[1].coeff == CB(1)
        @test res.poly.terms[1].coeff isa CB
    end
end

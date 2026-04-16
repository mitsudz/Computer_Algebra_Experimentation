using Test
using DynamicPolynomials

# Ensure the source file is included correctly
include("../multivar_poly_rep.jl")

@testset "FastPoly GrLex Representation Tests" begin
    @polyvar x1 x2 x3 x4 x5
    vars = [x1, x2, x3, x4, x5]
    
    @testset "Bit Packing and Ordering" begin
        m1 = pack_exponents([2, 0, 0, 0, 0]) 
        m2 = pack_exponents([1, 1, 0, 0, 0]) 
        
        @test m1.bits > m2.bits
        
        m3 = pack_exponents([0, 0, 0, 0, 3]) 
        m4 = pack_exponents([2, 0, 0, 0, 0]) 
        @test m3.bits > m4.bits
    end

    @testset "LCM Functionality" begin
        m_a = pack_exponents([2, 1, 0, 0, 0]) 
        m_b = pack_exponents([1, 2, 0, 0, 0]) 
        
        l = lcm(m_a, m_b, 5)
        unpacked = unpack_exponents(l, 5)
        
        @test unpacked == [2, 2, 0, 0, 0]
        @test (l.bits >> 112) == 4 
    end

    @testset "Divisibility" begin
        m_target = pack_exponents([2, 2, 0, 0, 0]) 
        m_div    = pack_exponents([1, 1, 0, 0, 0]) 
        m_notdiv = pack_exponents([3, 0, 0, 0, 0]) 
        
        @test divides(m_div, m_target) == true
        @test divides(m_notdiv, m_target) == false
    end

    @testset "Arithmetic (Addition & Multiplication)" begin
        dp_p1 = 2*x1^2 + 3*x2
        fp_p1 = from_dynamic(dp_p1, vars)
        
        t = Term(5, pack_exponents([1, 0, 0, 0, 0]))
        res = t * fp_p1
        
        dp_res = to_dynamic(res, vars)
        @test dp_res == 10*x1^3 + 15*x1*x2
        
        dp_p2 = x1^3 + x2
        fp_p2 = from_dynamic(dp_p2, vars)
        
        sum_fp = fp_p1 + fp_p2
        @test to_dynamic(sum_fp, vars) == (x1^3 + 2*x1^2 + 4*x2)
    end

    @testset "Leading Term Getters" begin
        p = from_dynamic(x1*x2 + x5^5, vars)
        @test leading_coefficient(p) == 1
        @test unpack_exponents(leading_monomial(p), 5) == [0, 0, 0, 0, 5]
    end
end

@testset "Advanced FastPoly Tests" begin
    @polyvar x1 x2 x3 x4 x5
    vars = [x1, x2, x3, x4, x5]

    @testset "Zero Polynomial & Cancellation" begin
        p1 = from_dynamic(x1 + x2, vars)
        p2 = from_dynamic(x1 + x2, vars)

        res = p1 + (Term(-1, pack_exponents([0,0,0,0,0])) * p2)
        @test iszero(res)

        @test length((p1 + (from_dynamic(-x1 - x2, vars))).terms) == 0
    end

    @testset "Rational{BigInt} Stress & Simplification" begin
        c1 = Rational{BigInt}(2^60, 3)
        c2 = Rational{BigInt}(1, 3)

        p = FastPoly([Term(c1, pack_exponents([1,0,0,0,0]))])
        t = Term(c2, pack_exponents([0,1,0,0,0]))

        res = t * p
        @test res.terms[1].coeff == Rational{BigInt}(2^60, 9)
        @test unpack_exponents(res.terms[1].mono, 5) == [1, 1, 0, 0, 0]
    end

    @testset "Cyclic-5 Generator Verification" begin
        poly_gen = x1 + x2 + x3 + x4 + x5
        fp_gen = from_dynamic(poly_gen, vars)

        @test length(fp_gen.terms) == 5
        @test unpack_exponents(fp_gen.terms[1].mono, 5) == [1, 0, 0, 0, 0]
        @test unpack_exponents(fp_gen.terms[5].mono, 5) == [0, 0, 0, 0, 1]

        poly_quad = x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1
        fp_quad = from_dynamic(poly_quad, vars)

        @test length(fp_quad.terms) == 5
        @test leading_coefficient(fp_quad) == 1
        @test unpack_exponents(leading_monomial(fp_quad), 5) == [1, 1, 0, 0, 0]
    end

    @testset "In-place Type Consistency" begin
        C = Rational{BigInt}
        p = FastPoly([Term(C(1, 2), pack_exponents([1,0,0,0,0]))])
        t = Term(2, pack_exponents([0,1,0,0,0])) 

        res = t * p
        @test res isa FastPoly{C}
        @test res.terms[1].coeff == 1 
    end
end

@testset "Recent Additions: Promotion and Monomial Mult" begin
    @polyvar x1 x2 x3 x4 x5
    vars = [x1, x2, x3, x4, x5]
    convert_poly = generate_converter(vars)
    
    @testset "Monomial to Poly Promotion" begin
        m = pack_exponents([1, 1, 0, 0, 0])
        
        fp_int = to_poly(m)
        @test fp_int isa FastPoly{Int}
        @test length(fp_int.terms) == 1
        @test fp_int.terms[1].coeff == 1
        
        fp_rat = to_poly(m, Rational{BigInt}(1, 2))
        @test fp_rat isa FastPoly{Rational{BigInt}}
        @test fp_rat.terms[1].coeff == 1//2
    end

    @testset "Monomial * Poly Multiplication" begin
        p = from_dynamic(x1 + x2, vars)
        m = pack_exponents([1, 0, 0, 0, 0])
        
        res = m * p
        @test to_dynamic(res, vars) == x1^2 + x1*x2
        
        res2 = p * m
        @test res2.terms[1].mono.bits == res.terms[1].mono.bits
    end

    @testset "Universal Converter" begin
        dp = x1^2 + x5
        fp = convert_poly(dp)
        @test fp isa FastPoly
        @test to_dynamic(fp, vars) == dp
        
        v_dp = [x1 + x2, x3 + x4]
        v_fp = convert_poly(v_dp)
        @test v_fp isa Vector{<:FastPoly}
        @test length(v_fp) == 2
        
        # This is where the error was. Just check the result directly.
        v_back = convert_poly(v_fp)
        @test v_back[1] == x1 + x2
        @test v_back isa Vector{<:DynamicPolynomials.AbstractPolynomial}
    end

    @testset "iszero helper" begin
        p_zero = FastPoly(Term{Int}[])
        p_val = to_poly(pack_exponents([0,0,0,0,0]))
        
        @test iszero(p_zero) == true
        @test iszero(p_val) == false
    end
end

@testset "FastPoly Scalar Multiplication" begin
    @polyvar x1 x2
    vars = [x1, x2]
    
    # Setup a basic poly: 2x1 + 4x2
    p_dynamic = (2//1)x1 + 4x2
    fp = from_dynamic(p_dynamic, vars)
    
    @testset "Integer Scaling" begin
        # (2x1 + 4x2) * 3 = 6x1 + 12x2
        res = fp * 3
        @test res isa FastPoly{Rational{Int}}
        @test to_dynamic(res, vars) == (6//1)x1 + 12x2
        
        # Commutativity: 3 * fp
        res_left = 3 * fp
        @test to_dynamic(res_left, vars) == (6//1)x1 + 12x2
    end
    
    #=
    @testset "Field Promotion (Int -> Rational)" begin
        # (2x1 + 4x2) * (1//2) = x1 + 2x2
        res_rat = fp * (1//2)
        
        # Verify the type evolved to handle the field
        @test res_rat isa FastPoly{Rational{Int64}}
        @test to_dynamic(res_rat, vars) == x1 + 2x2
        @test res_rat.terms[1].coeff == 1//1
        @test res_rat.terms[2].coeff == 2//1
    end
    =#

    @testset "Zero Multiplication" begin
        # Scaling by 0 should technically return an empty term vector
        res_zero = fp * 0
        @test iszero(res_zero)
        @test length(res_zero.terms) == 0
    end

    #=
    @testset "BigInt Promotion" begin
        # Ensure it plays nice with your thesis requirement for BigInts
        fp_big = fp * big(1)
        @test fp_big isa FastPoly{BigInt}
        
        res_big_rat = fp_big * (1//3)
        @test res_big_rat isa FastPoly{Rational{BigInt}}
        @test res_big_rat.terms[1].coeff == 2//3
    end
    =#
end

@testset "Identity Monomial Tests" begin
    # CRITICAL: Define polyvars inside the testset scope 
    # to ensure from_dynamic can resolve them.
    @polyvar x1 x2 x3 x4 x5
    vars = [x1, x2, x3, x4, x5]
    num_vars = length(vars)
    id = identity_monomial()

    @testset "Basic Properties" begin
        @test id.bits == 0
        @test (id.bits >> 112) == 0
        
        exps = unpack_exponents(id, num_vars)
        @test all(iszero, exps)
        @test length(exps) == num_vars
    end

    @testset "Multiplication Invariance" begin
        # Create a specific monomial: x1^2 * x3^1
        m1 = pack_exponents([2, 0, 1, 0, 0])
        
        # Test bit addition directly
        res_bits = m1.bits + id.bits
        @test res_bits == m1.bits
        
        # Test FastPoly * GrLexMonomial (identity)
        # Using the vars defined above
        p = from_dynamic(x1^2 + x2, vars)
        res_p = p * id
        
        @test res_p.terms[1].mono.bits == p.terms[1].mono.bits
        @test length(res_p.terms) == length(p.terms)
    end

    @testset "Divisibility" begin
        m1 = pack_exponents([1, 1, 1, 1, 1])
        @test divides(id, m1) == true
        @test divides(m1, id) == false
        @test divides(id, id) == true
    end
end

@testset "Total Degree Extraction" begin
    # 1. Identity Monomial (Degree 0)
    @test total_degree(identity_monomial()) == 0

    # 2. Specific Monomial: x1^3 * x2^2 (Total degree 5)
    # Bits layout: [5][3][2][0][0][0]...
    m1 = pack_exponents([3, 2, 0, 0, 0])
    @test total_degree(m1) == 5

    # 3. High degree check (testing 16-bit boundary)
    # Suppose a monomial has total degree 500
    exps = zeros(Int, 5)
    exps[1] = 500
    m_high = pack_exponents(exps)
    @test total_degree(m_high) == 500

    # 4. Consistency with addition
    # In GrLex bits addition: deg(m1 * m2) = deg(m1) + deg(m2)
    m2 = pack_exponents([1, 1, 1, 1, 1]) # Degree 5
    m_combined = GrLexMonomial(m1.bits + m2.bits)
    @test total_degree(m_combined) == 10
end

@testset "FastPoly Total Degree Tests" begin
    @polyvar x1 x2 x3
    vars = [x1, x2, x3]
    C = Rational{Int64}

    @testset "Basic Degrees" begin
        # 1. Zero polynomial
        # We use an explicit typed empty vector to avoid inference issues
        empty_p = FastPoly(Term{C}[])
        @test total_degree(empty_p) == 0

        # 2. Constant polynomial: f = 5 (Degree 0)
        # We multiply by one(x1) to ensure from_dynamic receives a
        # DynamicPolynomials object rather than a raw Rational.
        const_p = from_dynamic(C(5) * one(x1), vars)
        @test total_degree(const_p) == 0

        # 3. Univariate: f = x1^10 (Degree 10)
        uni_p = from_dynamic(x1^10, vars)
        @test total_degree(uni_p) == 10
    end

    @testset "Leading Term Dominance" begin
        # f = x1^2*x2 + x1*x2 + x3^2
        # Leading term is x1^2*x2 (deg 3), second is x1*x2 (deg 2), third is x3^2 (deg 2).
        # GrLex order ensures x1^2*x2 is at index 1.
        poly = from_dynamic(C(1)*x1^2*x2 + C(1)*x1*x2 + C(1)*x3^2, vars)

        @test total_degree(poly) == 3
        # Ensure it didn't mistakenly return the degree of the trailing terms
        @test total_degree(poly) != 2
    end

    @testset "After Multiplication" begin
        # Test that total_degree correctly reflects changes after monomial multiplication
        p = from_dynamic(x1 + x2, vars) # Degree 1
        m = pack_exponents([2, 1, 0])   # Degree 3 (x1^2 * x2)

        # res should be x1^3*x2 + x1^2*x2^2 (Degree 4)
        res = p * m
        @test total_degree(res) == 4
    end
end

@testset "GrLexMonomial Equality Tests" begin
    num_vars = 5
    
    @testset "Basic Equality" begin
        m1 = pack_exponents([1, 2, 0, 0, 0])
        m2 = pack_exponents([1, 2, 0, 0, 0])
        m3 = pack_exponents([2, 1, 0, 0, 0]) # Same total degree, different exponents
        
        @test m1 == m2
        @test m1 != m3
        @test m2 != m3
    end

    @testset "Identity and Constants" begin
        id1 = identity_monomial()
        id2 = GrLexMonomial(UInt128(0))
        m_const = pack_exponents([0, 0, 0, 0, 0])
        
        @test id1 == id2
        @test id1 == m_const
    end

    @testset "Edge Cases" begin
        # Test equality with different variables having the same exponent
        # x1^1 vs x2^1
        mx1 = pack_exponents([1, 0, 0, 0, 0])
        mx2 = pack_exponents([0, 1, 0, 0, 0])
        
        @test mx1 != mx2
        
        # Test maximum exponent boundary
        m_max = pack_exponents([MAX_EXP, 0, 0, 0, 0])
        m_max_copy = pack_exponents([MAX_EXP, 0, 0, 0, 0])
        @test m_max == m_max_copy
    end
end

@testset "Critical Pair Logic" begin
    @polyvar x1 x2 x3
    vars = [x1, x2, x3]
    num_vars = 3
    C = Rational{Int64}

    @testset "Simple Coprime Monomials" begin
        # f = x1^2, g = x2^2
        # LCM = x1^2 * x2^2. u = x2^2, v = x1^2
        f = from_dynamic(C(1)x1^2, vars)
        g = from_dynamic(C(1)x2^2, vars)
        
        u, v = critical_pair(f, g, num_vars)
        
        @test unpack_exponents(u, num_vars) == [0, 2, 0]
        @test unpack_exponents(v, num_vars) == [2, 0, 0]
    end

    @testset "Overlapping Monomials" begin
        # f = x1^2 * x2, g = x1 * x2^2
        # LCM = x1^2 * x2^2. u = x2, v = x1
        f = from_dynamic(C(1)x1^2 * x2, vars)
        g = from_dynamic(C(1)x1 * x2^2, vars)
        
        u, v = critical_pair(f, g, num_vars)
        
        @test unpack_exponents(u, num_vars) == [0, 1, 0]
        @test unpack_exponents(v, num_vars) == [1, 0, 0]
    end

    @testset "Total Degree Check" begin
        f = from_dynamic(x1^3, vars)
        g = from_dynamic(x1*x2, vars)
        u, v = critical_pair(f, g, num_vars)
        
        # LCM is x1^3 * x2 (deg 4)
        # u = x2 (deg 1), v = x1^2 (deg 2)
        @test total_degree(u) == 1
        @test total_degree(v) == 2
    end
end

@testset "GrLexMonomial Multiplication" begin
    num_vars = 3
    
    @testset "Standard Multiplication" begin
        # m1 = x1^2 * x3^1 (Total Deg: 3)
        m1 = pack_exponents([2, 0, 1])
        # m2 = x1^1 * x2^4 (Total Deg: 5)
        m2 = pack_exponents([1, 4, 0])
        
        res = m1 * m2
        exps = unpack_exponents(res, num_vars)
        
        @test exps == [3, 4, 1]
        @test total_degree(res) == 8
    end

    @testset "Identity Property" begin
        m1 = pack_exponents([1, 5, 2])
        id = identity_monomial()
        
        @test m1 * id == m1
        @test id * m1 == m1
    end

    @testset "Bit Boundary Integrity" begin
        # Test that adding exponents doesn't "bleed" into the next slot
        # unless there is an actual overflow (> 65535)
        m1 = pack_exponents([100, 100, 100])
        m2 = pack_exponents([200, 200, 200])
        
        res = m1 * m2
        @test unpack_exponents(res, num_vars) == [300, 300, 300]
        @test total_degree(res) == 900
    end
end

@testset "GrLexMonomial Ordering Tests" begin
    @testset "Degree Precedence" begin
        # x1^2 (Deg 2) vs x1^1 * x2^2 (Deg 3)
        # Even though x1^2 has a higher exponent in the first slot, 
        # GrLex says Deg 3 > Deg 2.
        m_deg2 = pack_exponents([2, 0, 0, 0, 0])
        m_deg3 = pack_exponents([1, 2, 0, 0, 0])
        
        @test m_deg3 > m_deg2
        @test m_deg2 < m_deg3
    end

    @testset "Lexicographical Tie-breaking" begin
        # Both Degree 3
        # m1 = x1^2 * x2^1
        # m2 = x1^2 * x3^1
        # m3 = x1^1 * x2^2
        m1 = pack_exponents([2, 1, 0, 0, 0])
        m2 = pack_exponents([2, 0, 1, 0, 0])
        m3 = pack_exponents([1, 2, 0, 0, 0])
        
        # x1^2 * x2^1 > x1^2 * x3^1 (x2 vs x3)
        @test m1 > m2
        # x1^2 * x2^1 > x1^1 * x2^2 (x1^2 vs x1^1)
        @test m1 > m3
    end

    @testset "Identity Comparisons" begin
        id = identity_monomial()
        m1 = pack_exponents([0, 0, 0, 0, 1]) # x5^1
        
        @test m1 > id
        @test id < m1
    end
end

using Test

@testset "FastPoly Zero Tests" begin
    # 1. Type-based creation
    z_int = zero(FastPoly{Int64})
    @test iszero(z_int)
    @test z_int isa FastPoly{Int64}
    @test length(z_int.terms) == 0

    z_rat = zero(FastPoly{Rational{BigInt}})
    @test iszero(z_rat)
    @test z_rat isa FastPoly{Rational{BigInt}}

    # 2. Instance-based creation
    p = to_poly(identity_monomial(), 5) # 5
    z_p = zero(p)
    @test iszero(z_p)
    @test z_p isa FastPoly{Int64}

    # 3. Arithmetic Identity (p + 0 = p)
    m = pack_exponents([1, 2, 0])
    p1 = to_poly(m, 10) # 10*x*y^2
    @test p1 + z_p == p1
    @test z_p + p1 == p1

    # 4. Constructor consistency
    @test FastPoly{Float64}() == zero(FastPoly{Float64})
end

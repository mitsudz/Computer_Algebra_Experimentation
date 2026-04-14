using Test
using DynamicPolynomials

include("../reduce_gb.jl")
# Assuming these are available in your environment
# include("multivar_poly_rep.jl") 

function rationalise(ps) 
    if ps isa Vector
        return [(big(1)//1) * p for p in ps]
    end
    return (big(1)//1) * ps
end

@testset "Groebner Helper Tests (Dynamic -> FastPoly)" begin
    @polyvar x y z
    num_vars = 3
    
    # Create the converter closure for our variable space
    converter = generate_converter([x, y, z])

    @testset "Basis Reduction (reduce_gb)" begin
        
        # --- Test 1: Redundant Elements ---
        # Ideal <x, y>. Set contains x^2 + y which is redundant because x | x^2.
        # DynamicPolynomials input -> Convert to FastPoly
        G1 = converter(rationalise([x, y, x^2 + y]))
        red1 = _reduce_gb(G1, num_vars)
        expected = converter(rationalise([y, x]))
        
        @test red1 == expected

        # --- Test 2: Tail Reduction ---
        # G = {x + y, y} -> Reduced should be {x, y}
        G2 = converter(rationalise([x + y, y]))
        red2 = _reduce_gb(copy(G2), num_vars)
        expected = converter(rationalise([y, x]))
    
        @test red2 == expected

        # --- Test 3: Normalization ---
        # 3xy + 6z -> xy + 2z
        G3 = converter(rationalise([3x*y + 6z]))
        red3 = _reduce_gb(copy(G3), num_vars)
        expected = converter(rationalise([x*y + 2z]))
        
        @test red3 == expected
        @test leading_coefficient(red3[1]) == 1

        # --- Test 4: Idempotency (G == reduce(reduce(G))) ---
        # Using a non-reduced but valid GB for the ideal <x^2-y, xy-x, y^2-y>
        G4_dyn = rationalise([x^2 - y, x*y - x, y^2 - y])
        G4 = converter(G4_dyn)
        
        first_pass = _reduce_gb(copy(G4), num_vars)
        second_pass = _reduce_gb(copy(first_pass), num_vars)
        
        @test first_pass == second_pass
        @test all(p -> leading_coefficient(p) == 1, second_pass)
    end
end

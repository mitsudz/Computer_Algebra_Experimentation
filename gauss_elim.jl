#=
Take a matrix over a commutative ring. We can apply do matrix multiplication on these. Using a naive approach
this takes $n^2$ multiplications for matrix times vector. For matrix times matrix this takes $n^3$ multiplications.

Also consider (AB)u where A,B are matrices and u is a vector. (AB)u is an n^3 computation, but A(Bu) is an n^2
computation. The order of operations here matters for cost - just like for the order of multiplying many matrices.

This is also relevant for computing homomorphisms in general.

Looking at Gauss Elim for matrices over fields. We can use this for computing determinants, solving systems and 
inverting matrices.

To get a determinant we only need to reduce to upper triangular form.
We reduce fully to solve a linear system - so we can always get the determinant whilst computing any system.
To solve for the inverse, we also just do the same reduction steps, but apply them to the identity matrix.
=#

using LinearAlgebra

"""
Returns the determinant and the inverse matrix.
"""
function gauss_elim(mat::Matrix{T})::Tuple{Rational, Union{Nothing, Matrix{Rational}}} where T <: Number
    # Assume matrix is square
    @assert size(M,1) == size(M, 2) "Matrix must be square"

    n = size(M, 1)
    det_ = Rational(1) # Rolling determinant
    inv_ = Matrix{Rational}(I, n, n)  # Rolling inverse
    println(typeof(inv_))

    # Naive gauss elim to get upper triangular
    for i = 1:n
        if M[i, i] == 0
            # Swap row to a non-zero one in position i
            swapped = false
            for j = i+1:n
                if M[j, i] != 0
                    inv_[i, :], inv_[j, :] = inv_[j, :], inv_[i, :]
                    M[i, :], M[j, :] = M[j, :], M[i, :]
                    det_ *= -1
                    swapped = true
                    break
                end
            end
            if !swapped # Zero column (from position i down - don't need to reduce this column)
                det_ = 0
                continue 
            end
        end
        
        # Get a 1 in position i, i
        c = M[i, i]
        det_ *= c
        inv_[i, :] .*= (1 // c)
        M[i, i:n] .*= (1 // c)
        
        # Start reducing from position i, i down
        for j = i+1:n
            inv_[j, :] .-= M[j, i] * inv_[i, :]  
            M[j, i:n] .-= M[j, i] * M[i, i:n]  
        end
    end
    
    # M is now upper triangular
    det_ == Rational(0) && return det_, nothing # No inverse
    for i = n:-1:2
        for j = i-1:-1:1
            inv_[j, :] .-= M[j, i] * inv_[i, :]
        end
    end

    return det_, inv_
end

M = [1 7 3; 2 5 3; 4 8 1]
#M = [1 7; 2 3]
# M = [1 0 0; 0 1 0; 0 0 1]
M = Rational.(M)

println("Matrix was:")
display(M)
println("In built determinant gives: $(det(M))")
println("In-built inverse was:")
display(inv(M))
det_, inv_ = gauss_elim(M)
println("Calculated determinant was: $det_")
println("Calculated inverse was:")
display(inv_)

#=
Time complexity analysis:

In the lecture it does it in general for an nxm matrix getting the determinant of the principal submatrix and taking 
the matrix to reduced row echelon form (and if m = n we get just the matrix, if m=n+1 we get solving a linear system
and if m=2n we get finding the inverse).
=# 


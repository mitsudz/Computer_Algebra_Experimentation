#=
Notes:

Consider $n$ by $n$ matrix over commutative ring $R$. We want to compute the determinant of $A$ (the matrix) 
as fast as possible. Basic definition involves cofactor expansion. Actually computing this way requires
$n!$ number of multiplications which is incredibly slow. 

Gentleman and Johnson provided a paper with a better method, they looked at the number of times you compute the 
determinant of each 2 by 2 submatrix. We compute the det of any 2 by 2 submatrix from the last two rows 
(or columns depending on how you expand), $n-2$ times. We can probably just do dynamic programming and massively 
speed this up. There are $n$ choose $2$ of these 2 by 2 submatrices from the last two rows. That pattern continues
and decreases as you go up in submatrix size. Their algorithm just involved going backwards from the submatrices 
first and not recomputing. Basically this is just a bottom up solution (ie dynamic programming - ASK PAUL - is this
commonly taken advantage of in computer algebra?). This isn't exactly dynamic programming - it goes from factorial to
exponential which still isn't great.

Of course, we could still try doing Gaussian elimination to handle this - although this had some complexity (that I
don't remember) that wasn't insignificant. The guassian elimination strategy is way better so far (polynomial time
so far). This was over integers specifically though, not over general rings - ASK PAUL - is there any reason we can't 
do this determinant calculation way faster or why we don't know whether we can? Over fields gaussian elimination works
because we have inverses - over rings we don't necessarily have inverses so we can't even perform gaussian elimination.

QUESTION - What is coefficient swell?

The technique of using rationals to solve integer matrix determinants is extending the ring of integers to it's field
of fractions then re-interpreting the result as an integer - this is possibly generalisable to other fields. This 
only works for integral domains because only then is the field of fractions well defined.

Berkowitz in 1984 showed how we can compute determinants in order $n^4$ operations. Kaltofen in 1982 used a 'baby step
giant step' method (WTF - ASK PAUL) to bring this down to order $n^{3.5}$.
=#

#= When $R$ is an integral domain...

One option as mentioned above is running Gaussian Elim in the fraction field. This requires computing gcd's however,
which is expensive.

I guess it makes sense that we can avoid the fractions since the answer has to be in the ring (since the definition
of determinant only uses ring operations). To avoid the fractions when cancelling e.g. the first column, we can 
compute the LCM of the entries pairwise and multiply each one to get that. Then we can cancel without ever having to
enter the fraction field. How do we compute the LCM without just doing GCD (isn't that what we're trying to avoid
ASK PAUL) - ie, LCM(a, b) = ab / GCD(a, b).
-> Apparently this is worse than just computing via the definition? Maybe because the coefficients inside the matrix
will likely explode? Yes - the size of the integers grows exponentially in $n$, but the determinant is bounded 
(Hadamard bound) by $\prod_{i=1}^n (\sum_{j=1}^n A_{ij}^2)^{1/2}$. If the numbers in the matrix are bounded by $B^m$
then the determinant is bounded by $n^{n/2} B^{nm}$. This approximately (using a large base) gives that the number of
digits of the determinant is controlled by $mn$ (with a lot of constants being omitted). This is much smaller than
exponential.
=#

#= Bareiss-Edmonds - IMPLEMENT THIS

This is essentially the same as doing Gaussian elimination, but we change the row operation slightly to be an exact
division and control the size of all the matrix entries.

We actually get the determinant of each k by k submatrix in the top left along the diagonal as we go through this 
procedure.

We multiply at each row by the current pivot, and divide by the previous pivot (initialised to be 1). The other numbers
also correspond to certain determinants.
=#

"""Bareiss Edmonds for integers (although it can be done for any integral domain)."""
function bareiss_edmond_det(M::Matrix{BigInt})::BigInt
    # Note - use views to be more efficient in Julia
    @assert size(M, 1) == size(M, 2) "Must be a square matrix"
    n = size(M, 1)

    A = copy(M)
    prev_pivot = 1 # Rolling pivot (previous determinant)
    for piv_idx in 1:n
        pivot = A[piv_idx, piv_idx]
        pivot == 0 && return 0 # Determinant of a submatrix is 0
        for row in piv_idx+1:n
            # Apply operation to each row, multiplying row by pivot value and 
            # subtracting pivot row times first entry before dividing by previous pivot (IMPORTANT)
            first_entry = A[row, piv_idx]
            A[row, piv_idx:n] .*= pivot
            A[row, piv_idx:n] .-= first_entry .* A[piv_idx, piv_idx:n] 
            A[row, piv_idx:n] .÷= prev_pivot
        end
        prev_pivot = pivot
    end

    display(A)
    return A[n, n]
end

M = [2 1 3; 3 -1 1; 5 3 1]
M = [213 2131 1231 9843; 2343 2234 1233 4136; 3324 3244 6442 2136; 2314 5643 6542 3456]
M = BigInt.(M)
#M = [2 1 3 4; 3 -1 1 1; 5 3 1 4; 3 6 7 2]
println(M)
println("Mine: $(bareiss_edmond_det(M))")
using LinearAlgebra
M = Rational{BigInt}.(M)
println("In built: $(det(M))")

# Is this accurate or is there a bug in mine?? - ASK PAUL, unsure exactly what is going on here

#= For polynomials

Apparently this isn't so good for polynomials. This is because at each step we multiply by the pivot, but
this pivot gets larger and larger the further along we go in the process. For integers, multiplying two entries
that are m and n digits gives us an integer of n + m digits. The degree of a polynomial is the same - ASK PAUL
about 1.23.45 in lec 6 - he says we may get 10000 terms for 100 degree polynomials multiplied together?
Maybe he means the number of intermediate terms could be up to 100^2?

Apparently the Gentlemen Johnson algorithm doesn't have blowup so will still be faster even though its algebraic 
complexity is exponential.

HE MEANS FOR MATRICES OVER MULTIPLE VARIABLES!
=#

#= Investigate symmetric matrices when this algorithm is run over this. 
=#

#= Testing whether the matrix of polynomials is singular?

Symmetric Toeplitz matrix is one with each diagonal from top left having the same values (and being symmetric).
Randomly selecting values in a matrix make it very unlikely to get a matrix with zero determinant.

Schwarz-Zippel Lemma (1978):
Let D be an integral domain, f a polynomial over D in n variables not equal to 0. Let S be a finite subset of D.
Suppose a chosen randomly from S. Then Prob(f(a) = 0) <= deg(f)/|S| where deg(f) is the total degree of f.

We apply this by letting the entries by variables in a matrix, then the polynomial f can represent the determinant, 
and we can bound it's total degree.

So to check whether a matrix is non-singular, we can use a probabilistic method, by simply picking a reasonably large
set S, and evaluating the polynomial that is the determinant k times (DON'T WE STILL HAVE TO GET THE DETERMINANT?)
and this gives us (deg(f)/|S|)^k probability of f being zero at all evaluation points.
 
ASK PAUL - VERY CONFUSED ABOUT HOW THIS WORKS

This helps us test if a matrix is non-singular for polynomials, since you can just test a few evaluation points,
and if you keep hitting roots the probability that happened due to random chance becomes exceedingly small.


APPLICATION - testing if a graph has a perfect matching (non-singularity of Tutte matrices or Edmonde matrices).
=#

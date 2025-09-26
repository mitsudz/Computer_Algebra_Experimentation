#= Adjugate Matrix

The adjugate matrix is the transpose of the cofactor matrix.
The cofactor matrix contains element c_{i, j} defined as 
(-1)^{i+j} \cdot det(A'_{i, j}) where A'_{i, j} is the n-1 x n-1
submatrix of A obtained by removing row i and column j (i.e. signed 
determined - very much connected to representation theory I would 
imagine). 

We have when the inverse is defined that adj(A) = A^{-1} \cdot det(A)
and the adjugate is defined in general (regardless of whether the inverse
exists) which makes it useful both in theory and practise. This property
also maintains the ring structure I think (ie clears field of fractions 
by mulitplying by det).

Another property that holds in general is that composition of 
adjugate and multiplication is almost commutative -
adj(AB) = adj(B) \cdot adj(A)

Pf:
To show adj(A) = A^{-1} \cdot det(A), first look at AX = I where X is the
inverse. Applying to the jth column of X, we get A \cdot x_j = e_j and we
apply Cramer's rule to the linear system over j's. det(A^(i)) has e_j
in column i, so it's determinant uses the determinant of the matrix A
removing column i and row j. Then we get the sign (-1)^{i+j} by going down
the row, and we get the determinant of the submatrix which gives
us the c_{j, i} cofactor. We then get that x_{i, j} = c_{j, i} / det(A) 
using Cramers rule. But c_{j, i} = adj(A)_{i, j} and we get the result.

Cramers Rule: Let Ax = b where A is square. Set A^(i) to be A with column
i replaced by b (in our case e_i). Then x_i (in our case X_{i,j}) is 
det(A^(i)) / det(A). 
=#

#= Characteristic Polynomial
Given a square matrix over a ring, the characteristic polynomial
c(y) := det(A - yI) which is a polynomial in R[y]. The eigenvalues
of A are the roots of c(y).

Berkowitz Algorithm: Computes c(y) using O(n^4) ring operations.
Note, c_0 = c(y=0) = det(A - 0I) = det(A) - this is a byproduct of the 
algorithm. The algorithm partitions the matrix first. Let Ar be the 
r x r principal submatrix and R, S be vectors in the ring of size R (R
along the bottom, S along the right),
with the bottom right element being A_{r+1, r+1}. The algorithm will 
recursively calculate char(Ar) changing r at each step and then combine
this to create char(A).

Thm 1: det(A) = det(Ar)A_{r+1, r+1} - R^T \cdot (adj(Ar) \cdot S). 

Thm 2: Let c_r(y) = det(Ar - yI) = sum_i c_i y^i
Then adj(Ar - yI) = - sum_{k=1}^r sum_{j=0}^{r-k} c_{k+j} Ar^j y^{k-1}

How do we get this, who knows? ASK PAUL

Berkowitz Cont:
Let A be an nxn matrix over R. Apply Thm1 to A - yI, with r = n-1
cn(y) = det(A- yI) = det(Ar - yIr)(A_{n, n} - y) - R^T \cdot (adj(Ar - yIr) \cdot S).
We now manipulate this expression. det(Ar - yIr) is cr(y) (recursive call). We apply
Thm2 to the adjugate of Ar - yI (where c_{i} is the coefficient of the char poly for 
cr(y)). \cdot is either matrix/vector multiplication or dot product where it makes sense.

So cn(y) = cr(y) (A_{n, n} - y) - R^T \cdot (-sum_{k=1..r} sum_{j=0..r-k} c_{k+j} Ar^j y^{k-1}) \cdot S
= cr(y) (Ann - y) + sum_{k=1..r} sum_{j=0..r-k} c_{k+j} (R^T \cdot Ar^j \cdot S) y^{k-1}

This is the entire algorithm, with r=n-1. The coeffiecients are obtained by computing cr(y), and the element
R^T \cdot Ar^j \cdot S goes to a scalar in the ring in terms of j so we denote them as Qj.

So, cn(y) = cr(y) (Ann - y) + sum_{k=1..r} sum_{j=0..r-k} c_{k+j} Qj y^{k-1}

We first compute Qj once, and compute the coefficients c_i once. The j's can range from 0 to r-1 (when
k=1).

Algorithm to compute Qj: (\cdot is dot product)
Qj = R^T \cdot (Ar^j S)

B <- Ir
Q0 <- R^T \cdot S // r mults
for j = 1..r-1 do
    B <- Ar \cdot B // = Aj^r, r^3 mults
    Qj <- R^T \cdot (B S) // r^2 + r mults

In total this is of order r^4 (which is n^4)

There is a faster method which avoids matrix matrix mults:
We do R^T ( Ar^{j-1} (Ar S))
Given, R, S, Ar ==>

Q0 <- R^T \cdot S // r
for j = 1...r-1 do
    S <- Ar S // = Ar^j S, r^2
    Qj <- R^T \cdot S // r

This is of order n^3.

Final Algorithm:
Input: A \in R^{nxn}, R a commutative ring. 
Output: c(y) := det(A - yI) \in R[y]

Berkowitz(A) {

if n = 1: output A11 - y

Let A be partitioned into Ar, R, S, Ann

Compute cr(y) = Berkowitz(Ar) = sum_{i=0..r} c_i y^i

Q0 <- R^T \cdot S
for i = 1...r-1 do {
    S <- Ar S
    Qj <- R^T \cdot S
}

output cr(y) (Ann - y) + sum_{k=1..r} sum_{j=0..r-k} c_{k+j} Qj y^{k-1}
}


Note - we do need the ring to have an identity so we can have identity matrices, so 
R needs to be a commutative unital ring.

For time complexity: M(n) = M(n-1) + O(n^3) = O(n^4)
M(n) = 1/4 n^4 - 1/3 n^3 + ...

Sparse Matrices: (Remark)
Gaussian Elim doesn't work well for sparse matrices since it can't recognise this.
Berkowitz does still work well in that case.
If A is sparse, the matrix multiplications don't have to be so slow if coded carefully.
We could do it linearly if the number of non-zero terms is linear in the dimension.
This for Berkowitz would speed it up by a factor of n!
=#


#= SUMMARY
To compute det(A) over R:

R is a field: Ordinary Gaussian Elimination (1/3 n^3 + ... mults in R)
R an integral domain: Bareiss/Edmonds does 2/3 n^3 + ... mults in R (and rational reconstruction)
    and 1/3 n^3 exact divisions in R
R is a commutative unital ring: Berkowitz does 1/4 n^4 + ... mults in R

Extension - Kaltofen's algorithm: order n^3.5 mults in R for R a commutative unital ring.
=#

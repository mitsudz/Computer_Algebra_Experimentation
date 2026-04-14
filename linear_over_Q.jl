#= Intro:
Want to solve Ax = b (e.g., A, b over Z and x over Q solving for x).

Cramer's rule:
Take Ax=b, A, b integer and x over Q. Replace ith column of A by b to get A^i. Then
x_i = y_i/D where y_i is the determinant of A^i and D = det(A).

How big can y_i and D be if all entries are bounded by B^m. Apply the Hadamard bound from here (see lec 7 notes).
=#

#= Bareiss/Edmonds/Dickson Algorithm
Computes y_i and D with no fractions, then creates x_i by fractions and computes gcd's there.
Computing these takes order n^3 ring operations (exact division included) of integers of size
order mn, where maximum input matrix entry size is B^m.

In terms of bit operations this is quite slow since these ring operations are all quadratic operations in 
the length of the entries which can go up the mn in size, so overall is order n^5m^2 in bit operations.

Essentially we triangularise the matrix using only integer operations. Uses something called forward and back 
elimination steps.
=#

#= SIDE NOTE - IMPORTANT FOR NEXT ALGORITHM!
Congruence Law:

***
ac = bc mod n and d = GCD(c, n) => a = b mod n/d
***
 
This in particular let's us modify things using the p-adic approach!
=#

#= Dixon/Monck Algorithm - DEFINITELY TRY TO IMPLEMENT THIS!

Improves cost (in bit operations) to n^3m^2 (2 orders of magnitude).

Solve Ax=b mod p^k with p large enough so x can be recovered using `rational number reconstruction`

Want to write x in the base p representation (p-adic rep).

Step 1: Construct A^{-1} mod p. If p | det(A) so A^{-1} doesn't exist, try another prime.
If p | det(A) for a lot of randomly chosen primes, it is likely that the matrix is singular.
Note, since we have a bound on the determinant, we could try enough primes such that the LCM of those
primes is larger than the bound proving it must be 0.

Step 2: Solve (X) b - Ax = 0 mod p^k for x for some k. Write x in base p representation.
So x = x_0 + x_1p + ... + x_{k-1}p^{k-1}
Solving (X) mod p gives b-Ax_0 = 0, so x_0 = A^{-1}b mod p.
Then solve mod p^2, which gives b - A(x_0 + x_1p) = 0 mod p^2
which gives => (b-Ax_0) - Ax_1p = 0 mod p^2
By the earlier part, (b-Ax_0) is divisible by p, so (apparently there are some rules around congruences where we
can divide like this:)
=> (b-Ax_0)/p - Ax_1 = 0 mod p
So x_1 = A^{-1} \cdot (b-Ax_0)/p
See lec notes for the next term...
Note, for the next term even though we divide by p^2 in theory, we don't want to let the powers
of p increase in practise. Thus, we reuse the previous division by p to split up the division
into only doing a division by p.

Step 2 Pseudocode:
e0 = n, x = 0^n
for k = 0,1,...
    xk = A^{-1} ek modp
    x += xk * p^k
    # Check if we are done
    if k+1 is a power of two:
        Let y = result of rational number reconstruction of x mod p^k
        If y didn't fail and b-Ay = 0 then outut y//
    ek+1 = (ek - A xk)/p # Exact division


Additional notes:
1) Only compute A^{-1} mod p once - this is order n^3
2) e{k+1} = (ek - A \cdot xk)/p is all exact division and integer multiplication, and is order n^2m (the most 
expensive step of the algorithm). Also xk is over mod p so is relatively `small` numbers.
3) We can do the check with b-Ay = 0 using integer arithmetic by computing the LCM=l of components of y and
checking if l*B - A(l*y) = 0.
=#

#= Rational Number Reconstruction 
Suppose n/d is in Q, and n,d in Z with gcd(n, d) = 1, d not 0. Suppose we have computed u = n/d mod m where
0 <= u < m and gcd(m, d) = 1 (otherwise existence fails). 

The problem is recovering n/d from u mod m...

We need to add an additional constraint to get uniqueness. We need m > 2|n|d (the 2 is for the +/-). This is
the smallest that can work.

We run Ext. Euc. Alg. with input m>u>=0.
We obtain si, ti, ri such that si*m + ti*u = ri for 0<=i<=N+1 where r{N+1} = 0

mod m we get ti*u = ri mod m
So u = ri/ti mod m if gcd(ti, m) = 1 (if we chose a prime this would work most times).
Cannot work for i = 0, since t0 = 0, also i < N+1 since t{N+1} = m = 0 mod m and so is not invertible.

So any 0 < i < N+1 could be a solution that we want!
Is n/d = ri/ti for some i - and how do we get such an i?

Theorem: n/d = ri/ti for some i if m > 2|n|d!

Note - the solution may be one of these ri/ti for a smaller m than this condition, but this condition
guarantees that we get the result we wanted.

Theorem (Euy, Davenport, Wang) 1982:
Let n, d in Z, d > 0, gcd(n, d) = 1. Let m in Z, m > 0, gcd(m, d) = 1.
Let u = n/d mod m with 0 <= u < m
Let N >= |n|, D >= d. Then,
(1) if m > 2ND, phi_m(n/d) is one to one (essentially the bound gives us enough numbers to work with). phi_m means
reduce mod m -> so there are enough integers mod m for fractions satisfying the bound. 
(2) If m > 2ND, then on input of m, u there exists a unique index i in the EEA s.t. ri/ti = n/d. Moreover, i is the
first index s.t. ri <= N.

Note - when running this algorithm we really need to know how big we want the numerator around denomenator to be.
If we have good bounds we can just get to a solution with a sufficient modulus, then run rational reconstruction.
However, if we don't have good bounds we'll just have to keep trying (as we do in the previous algorithm).
We don't just use a large modulus immediately, since it's possible that the solution is relatively small.

For solving Ax = b:
Wang says just try N = D = floor(m^{1/2}/2) and try rational reconstruction. This may NOT succeed, since we may get
that the first ri < N doesn't have ti < D which is a failure (which means we need a bigger modulus).

One thing to notice, is for the solution i, the quotient is relatively large and every other quotient is much smaller.
Monagan proved (2003): |n|*d*q_{k+1} > m/3 -> this shows if n and d are small relative to the modulus, then the 
quotient must be large (and vice versa).

This means that you can just select the row with the largest quotient to get the solution. It also somehow means
that you can use maximal quotient RR to use smaller moduli and get the reconstruction faster.
TODO - Investigate this further! The lecture is a very brief introduction to it!
=#


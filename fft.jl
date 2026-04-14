#= The first version of FFT displayed is the Cooley-Tukey Implementation

First idea: how to evaluate and interpolate polynomials in one variable?
Horner's and Newton's methods give us O(n^2). FFT gives us O(nlogn) for 
both sides fast. Consider univariate polynomial over a ring.

Let a = [a_0,...,a_{n-1}] be a polynomial in x, and y1,..,yn points to 
interpolate at, and [a(y1),..., a(yn)] another vector. 
Idea of FFT is to first split up a(x) into even and odd powers:
a(x) = a_0 + a_1 x + ... + a_{n-1} x^{n-1} = 
(a_0 + a_2 x^2 + ... + a_{n-2}x^{n-2}) + x(a_1 + a_3 x^2 + ... + a_{n-1} x^{n-2})
assuming n is even.
Then set the first part to be b(x), and the second to be c(x) so we get:

a(x) = b(x^2) + x*c(x^2)

This gives us two polynomials with half the terms, and so evaluating via Horner's method
takes half the time, and by evaluating at +/- y we get the same b(y^2), c(y^2) and so
only need to evaluate once to get a(y) and a(-y) effectively cutting out half the work. We
only really had to do one extra multiplication.

We can extend this idea by utilising roots of unity in our field (in most applied fields,
we use complex numbers to do this). We have the usual facts about roots of unity from ring 
theory (e.g, sum_{i=0..n-1} w^i = 0 where w is a primitive nth root of unity). This needs
to be proven in general rings (or is it fields?), but some results hold beyond just complex numbers.

DEF: Let w be a primitive nth root of unity in a field F. Given a polynomial of degree <= n-1, 
let B := [a(w^0=1),...,a(w^{n-1})] is the discrete Fourier Transform of a(x).

How fast can we compute B? Just doing Horner's would give us quadratic time.

Algorithm: Discrete Fast Fourier Transform (DFFT):
Inputs: (over a field F)
    n = 2^k
    A = [a_0,...,a_{n-1}], a(x) = sum_i a_i x^i
    w an nth root of unity

Output: [a(1), a(w),..., a(w^{n-1})]

if n = 1 {
    // A = [a_0]
    // a(x) = a_0
    // a(1) = a_0
    return A
}

b <- [a_0, a_2,...,a_{n-2}] // b(x) = a_0 + a_2 x + ... + a_{n-2} x^{(n-2)/2}
c <- [a_1, a_3,...,a_{n-1}] // c(x) = a_1 + a_3 x + ... + a_{n-1} x^{(n-1)/2}
// a(x) = b(x^2) + x c(x^2)

B <- DFFT(n/2, b, w^2) // B = [b(1), b(w^2), b(w^4),...,b(w^{n-2})]
C <- DFFT(n/2, c, w^2) // C = [c(1), c(w^2),..., c(w^{n-2})]

y <- 1 // Recording powers of w
for i = 0,1,..,(n/2)-1 do
    T <- y C[i] // y = w^i
    A[i] <- B[i] + T // A[i] = b(w^{2i}) + w^i c(w^{2i}) = a(w^i)
    A[i + n/2] <- B[i] - T
    y <- w y
end for // Cost of each loop in the field is 2 multiplications and additions

return [A_0, A_1,...,A_{n-1}]

COST: 
n=1 -> T(1) = 0
n>1 -> T(n) = 2T(n/2) + 1 + n // extra 1 for w^2, and rest is for loop

Solving this recurrence gives T(n) = nlog_2 n + n-1 which is nlogn time.

Optimisation: Precompute 1,w,..,w^{n/2 - 1} using n/2 - 1 mults in F
Set W = [1,w,..,w^{n/2 - 1}, 1, w^2, w^4,..., w^{n-2}, 1, w^4,...,w^{n-4},...,1,0]
These are the values used in the recursive calls at each level, 
(with the padded 0 this has length n) - by simply passing sections of this array
down we never need to compute powers of w.

We then only have one multiplication in the loop, and the actual time complexity changes to:
n>1 -> T(n) = 2T(n/2) + n/2 // Also don't need to compute w^2
with a size n/2 overhead to compute the powers intially (and size n to actually store them)

We then get T(n) = (nlog_2n)/2
=#

#= Inverse Transform:
Let a(x) = sum_{i=0,..,n-1} a_i x^i, n = 2^k and w a primitive nth root of unity in F.
Let A = [a_0,..,a_{n-1}]
Let B = [a(1), a(w),..., a(w^{n-1})]

Observe the (Van DeMonde) matrix V_w where each row is the 
powers of omega in increments of i for row i (indexing from 0),
multiplied by A^T gives the vector B^T.

For the inverse transform we aim to compute A from B. One way to compute B
would be to multiply V_w (van demonde matrix of w) by A. This takes n^2 mults.
One way to compute A (interpolate a(x)) is to solve V_w A = B for A.
Another way to compute A is to compute (V_w)^{-1} then A = (V_w)^{-1} B.
These two cost n^3 mults in F. For V_w, in position i, j (indexed from 1) we have
w^{(i-1)(j-1)} 

We can use linear algebra tools:

Lemma: Let w be a primitive nth root of unity. Then w^{-1} = w^{n-1} and 
w^{n-1} is also a primitive nth root of unity. Pf - easy ring theory stuff

Lemma: V_w \cdot V_{w^{-1}} = nI. So the inverse of V_w is V_{w^{-1}} / n
Partial Proof in lecture notes.

Thus, we don't actually need to compute the inverse, we can simply infer it.
So, A = (V_w)^{-1} B = 1/n \cdot V_{w^{-1}} B = 1/n DFFT(n, B, w^{-1})
since we can just pretend B is a polynomial!

So we can apply DFFT to both evaluate and interpolate (ie forward and backwards
transform).
=#

#= FFT Multiplication
Input: a,b in F[x], F a field
Output: c = a * b

Let n = 2^k be the first power of 2 > deg(c) = deg(a) + deg(b)
Find w in F primitive nth root of unity.
Let da = deg(a), db = deg(b)

A <- [a_0,...,a_{da},0,...,0]
B <- [b_0,...,b_{db},0,...,0]
such that A, B have length n.

A <- DFFT(n, A, w) // [a(1), a(w),...,a(w^{n-1})]
B <- DFFT(n, B, w) // [b(1), b(w),...,b(w^{n-1})]

// c(x) = a(x) * b(x)
C <- [A[1]*B[1],...,A[n]*B[n]]

// Now we can interpolate C
C <- 1/n * DFFT(n, C, w^{-1})
Output sum_{i=0,..,n-1} C[i] x^i

COST: 3 DFFT's of size n = 2^k > deg c and n multiplications
to compute C (as the product of A, B elementwise) and n for scaling
after the final DFFT.

So overall cost is 3/2 nlog_2n + 2n + n (last n is for computing powers of w and w^{-1})
so overall order is nlogn in F.

NOTE: We still haven't done anything about computing primitive nth roots of unity.
In complex numbers this is easy, in other fields not necessarily.

In Z_p w exists iff 2^k = n | p-1 (these are special primes - they have a name!) 
Assuming this:
1) Let e be a primitive element, ie e^{p-1} = 1 and e^k != 1 for 0<k<p
2) e^{p-1} = 1 => e^{nq} = 1 (where p-1 = nq for some q) so (e^q)^n = 1
and e^q is an nth root of unity and is primitive.

Getting a primitive element is another problem entirely.

=#

#= Other Optimisations & Implementations of FFT:
The current construction uses temporary arrays for B/C which overall takes up nlogn space.
We can modify our construction to only use O(n) space
We can provide to DFFT a temporary array of size n which will be enough
for all temporary storage for the DFFT to run in place.
We can then use some pointer magic to say where the array starts and let the temporary
array be T = B | C (concatenation).
Note we can do the FFT recursive calls in parallel.

Different version of FFT: (Section 8.2 of Modern Computer Algebra)
Let a(x) = sum_{i=0,..,n/2-1} a_i x^i + sum_{i=n/2,...,n-1} a_i x^i
(1) a ÷ (x^{n/2} - 1) => a(x) = q_0(x) (x^{n/2} - 1) + r_0(x), r_0(x) = a(x^{n/2} = 1)
r_0(x) = (a_0 + a_{n/2}) + (a_1 + a_{n/2 + 1}) x + ... + (a_{n/2 - 1} + a_{n-1}) x^{n/2 - 1)

(2) a ÷ (x^{n/2} + 1) => Similar to before, a(x) = q_1(x) (x^{n/2) + 1) + r_1(x),
r_1(x) = a(x^{n/2} = -1) which gives:
r_1(x) = (a_0 - a_{n/2}) + ... + (a_{n/2 - 1} - a_{n-1}) x^{n/2 - 1}

Then considering a(w^{2i}) and a(w^{2i+1}) offers some further insight, since
a(w^{2i}) = r_0(w^{2i})
a(w^{2i+1}) = r_1(w^{2i+1})
Setting r_1*(x) = r_1(w x) we get polynomials r_0 and r_1* of half the size and we need to evalutate
them at even powers of w so we can recursives call our FFT procedure as such.

See lec notes for how this is turned into pseudocode.
This is a genuinely different version of FFT, but it has the exact same cost in terms
of field operations.


Looking at what actually happens to the arrays in the C implementations, it appears we make certain
permutations of the arrays. When you do this, and you look at the permutation, you actually just reverse
the binary numbers:
Ie, if we started with [a0, a1,...,a7], we would end up with (in the original version of FFT):
[a0, a4, a2, a6, a1,a5, a3, a7] and if you write the subscript in binary, this is the same as 
writing 0-7 in binary and then reversing those representations! We call this the bit reversed
permutation and we can compute that in linear (in the number of bits) time!

This permutation is an involution!!! (Surely this must be related to how FFT works as being an involution
linear operator - is this somehow the FFT in some way in group theory land?). The number of moves on
this array is nlogn - this is expensive, can we eliminate this?

Looking at the second version of the FFT, the permutation happens after doing all the recursive calls (whereas
in the first version the permutation happens before the recursive call and the work is actually done).
In this version (see the lec notes) at the base of the recursion if there are only 2 elements in A, no actual
swapping occurs.

We could optimise the second version, since none of the computation occurs after the recursive call,
we can simply do the rearranging at the very end after doing all the recursive calls and computation
thus doing the reordering in linear time? (The order isn't exactly the same as before though).

I THINK both FFT versions do the same bit reversed permutation in some way, since they both rearrange
the first/last half with the even/odd! In fact they do the exact the same permutation! So in the second
version we don't need to permute at all!


For multiplication we can take advantage of the order of permutations to do polynomial multiplication even
faster. Given polynomials a, b - if we apply the second FFT, we do computation then permute to get A, B. 
Then after multiplying to get C = A*B elementwise we permute again then do computation if we apply
the first FFT. But in this case, the permutations would cancel since it is an involution, so we can
just not permute at all! So we apply the second FFT to a, b without permuting, then multiply A, B elementwise
to get C and then apply the first FFT without permuting to get c! This cuts out all the permutation work!
 
In fact, since we don't need to do the permutations, we don't even need the extra temporary array!

We now have 2 representations for our polynomial, a coefficient representation and a points representation.
The points representation makes it trivial to do multiplication, and we can convert between the two in nlogn
time. In the coefficient representation differentiation is easy but multiplication is not.

In addition, since the Fourier transform is a linear transformation, we can write
Fw(xA + yB) = xFw(A) + yFw(B) where w is a primitive nth root of unity greater than the degree of A and B,
and A, B are polynomials in coefficient or points representation and x, y are scalars in F our field.
We can take advantage of this.

E.g, suppose we want to multiply a matrix of polynomials by a vector of polynomials, we can do all the forward
transformations first, then apply the inverse transforms - see lec notes for details, but this saves a number
of Fourier transforms.


NOTE - THIS MATRIX TIMES POLYNOMIAL THING IS WHAT HAPPENS INSIDE THE FAST GCD ALGORITHM - TRY THIS ONE OUT!!!
=#

#= Fast Division
Let a, b in F[x] with deg(a) = m >= n = deg(b).

Let a = bq + r, deg(r) < deg(b). (This is a ÷ b).
The classical ÷ algorithm does (m-n+1)-n mults in F.
If m = 2n => (n+1)n is n^2 in n. 
We want an algorithm that is nlogn

(1) Compute q
(2) Compute r := a - b*q with fast multiplication

Define a^r = a_0 x^m + a_1 x^{m-1} + ... _ a_m the reciprocal polynomial.
Idea (1), compute q^r = a^r/b^r as a truncated power series to order m-n+1 which is the degree of q + 1
Then q = (q^r)^r
Note, division of power series is almost the same as normal division, except we start at the smallest power term
    and work from there up. This is genuinely a different idea.

Example in lec notes.
This long division to get up to the degree m-n+1 is still quadratic however. This is essentially the 
same as the normal division algorithm.

Idea (2) We want to compute a^r/b^r. Instead, compute 1/b^r to degree m-n+1 then 
compute q^r = 1/b^r * a^r to order m-n+1 using a second fast multiplication in F[x].

We can get the reciprocal faster: (Note just doing power series division again would still be n^2).
Theorem (Newton Iteration for Poly): 
Let R be a comm ring (F in our case). Let f be in R[x]. f = f0 + f1x + .. with f0^{-1} in R.
Let y0 = f0^{-1} and y_i = 2y_{i-1} - f*y_{i-1}^2 mod x^{2^i} for i > 0.
Then f*y_i = 1 mod x^{2^i} for i >= 0. (ie, yi = f^{-1} up to order x^{2^i}).

We want to compute y = 1/b (ignoring reciprocal notation) so we use f(y) = b - 1/y and attempt to solve
f(y) = 0 using Newton iteration. f'(y) = 1/y^2. So y_{k+1} = y_k - (b-1/y_k)/(1/y_k^2) = y_k - by_k^2 + y_k
= 2y_k - by_k^2 (which is where the the formula in the theorem comes from). In particular, we can compute
this without doing any divisions. This does require 2 more FFT multiplications, and then we need to truncate our
result to the right degree. We will start with y_0 = 1/b_0.

This approach doubles the order of our estimation at each step.
Let cost of multiplying 2 poly of deg <= n be m(n). Let I(n) be the number of arithmetic ops for computing 1/b^r mod
x^n using Newton iteration.
I(1) = 1. I(n) <= I(n/2) + m(n/2) + m(n) + cn 
// recursive call, squaring prev term, multiplying by b^r, scalar operations

Exercise: I(n) < 3m(n) + c'n
When we use FFT, this gives order nlogn
Let D(n) be cost of dividing a(x) ÷ b(x) to get q, r.
Then, D(n) = I(n) + m(n) + m(n) + cn // reciprocal, quotient, remainder
So D(n) = 5m(n) + c'n

So we can divide two polynomials in about as much time as multiplying 5 polynomials.

It is possible using the 'middle product' of Zimmerman to reduce I(n) < 2m(n), so D(n) <~ 4m(n). TRY THIS!!!

=#

#= Fast Multi-Point Evaluation
Assuming fast multiplication and division.
Let f in F[x], a1,..,an where n=2^k and deg(f) < n.
How fast can we compute all f(ai)?
Horner's method does <= n(n-1) (order n^2) mults and additions.

Since these are arbitrary points we can't just do FFT directly.
Consider the product tree Pi_n, for n=8 (example). This tree starts with the individual terms x-a_i,
and multiplies pairs of terms at each level. So level 1 has 8 terms of degree 1, level 2 has 4 terms of degree 2
and so on.

We write M_i to be the product in the nodes of this product tree, in BFS order (ie M_1 is the whole product
prod_{i=1,..,n} (x-a_i)).

Assuming we have computed the whole product tree, we can compute
r2 = f mod m2, r3 = f mod m3
r4 = r2 mod m4, r5 = r2 mod m5, r6 = r3 mod m6, r7 = r3 mod m7
r8 = r4 mod (x-a1), r9 = r4 mod (x-a2), r10 = r5 mod (x-a3) etc...

It turns out, r8 = f(a1), r9 = f(a2) etc...

We can do all our mults and divisions via FFT to get a good algorithm.

Let T(n) be the no of arithmetic ops of all the divisions. Suppose D(n) <= 4m(n) (the improved div alg - 
middle prod).
Then, T(n) = 2D(n/2) + 4D(n/4) + ... + nD(1) <= 4(2m(n/2) + 4m(n/4) +...+ nm(1))
Assuming multiplication is super linear (ie, at least as bad as linear time) then 
T(n) < 4logn m(n) // each term goes to m(n) by doing m(n) >= 2 m(n/2) and we have logn terms.

So T(n) is order m(n)*logn.
Note, if deg(f) >= n then we need to compute r1 = f mod m1 = prod_{i=1,..,n} (x-ai) which takes one extra mult
and ÷ (doesn't change overall complexity).

Lemma: Let f,g,h in F[x] where g|h. 
Then f mod g = (f mod h) mod g
Proof in Lec notes (pretty simply proof)

This lemma shows correctness for the product tree by being able to mod out by other larger factors first.

The cost of computing the product tree is as follows:
Assume we use FFT for all mult in F[x] in Pi_n. Then we have a problem to do with the degree of the polynomials
vs size of input to FFT!

If deg(fg) < n = 2^k we can multiply using 3 FFT's of size n. But in Pi_n, deg(fg) = n = 2^k. We need 3 FFT's of
size 2n. Exercise: Show how to do this using 3 FFT's of size n. (There will be a collisison somewhere and 
we just need to resolve that).

Let T(n) be no. of arithmetic ops in F needed to compute Pi_n. Let m(n) be for multiplication.
T(n) < n/2m(1) + n/4 m(2) + ... + 2 m(n/4) // (We don't need to compute the final multiplication in general)
Assuming again 2m(n/2) > m(n) we get that
T(n) < ... < M(n/2) + ... + M(n/2) = M(n/2) (logn - 1) < M(n) (logn - 1)/2
This is of order m(n)logn


So the total cost to compute f(ai) for all i=1,..,n is 2 * order m(n)logn which is m(n)logn
Other things in order m(n)logn include gcd, computing cube and square roots, computing inverse of a power series.
So multiplication is the bottleneck for a lot of fast algorithms.

=#


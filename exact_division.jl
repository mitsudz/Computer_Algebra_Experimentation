#=
Exact Division in F[x]:

Given f, g \in F[x], we want q \in F[x] s.t., f = q * g.
Define reciprocal of polynomial by 
reciprocal: p -> p_r
where p_r(x) = x^deg(p) * p(1/x)
Note: this is the same as reversing all the degrees of terms (ie, if deg(p) = 5, then the x^3 term becomes an x^2 term)

Then, f_r(x) = x^deg(f) * f(1/x) = x^{deg(g) + deg(f)} * g(1/x) * q(1/x) = g_r(x) * q_r(x)

Is reciprocal a ring endomorphism? Does it respect much of the structure
beyong multiplicative? Clearly it is it's own inverse.

To find q_r(x) we can use Newton iteration since we know the degree of q.
Newton iteration involves successively finding a polynomial p mod x^k for
successively larger k.
We want to find the inverse of g_r mod sufficiently large x^k, then
just multiply f_r by g_r^-1.

h_0 = 1 / const(g_r) where const is the constant term (cannot be 0 unless g was 0).
h_{k+1} =h_k * (2 − g_r * h_k) mod x^{2^k} // NOTE - mod x^{2^k} means discard terms of too high degree so we simply don't even need to compute them.

We just need to do this enough times to get k large enough.
For exact division this works perfectly with no errors in coefficients that
we find.

We then after computing f_r * (1/g_r) (e.g., using FFT or some other 
algorithm) we then just reverse the reciprocal q_r by taking q = reciprocal(q_r).

This gives us the time complexity being 2 reciprocals, one large multiplication
and several smaller multiplications to get the inverse of g_r. Apparently
this comes out to just O(M(n)) where M(n) is the cost of multiplying
two polynomials of size n.

Note - it seems the only reason we need to take the reciprocal is to ensure
the constant term is invertible (non-zero in F). If we already knew that
the constant term of g was non-zero we could just apply Newton iteration
to 1/g directly with h_0 = 1 / const(g).
=#


# TODO - IMPLEMENT EXACT DIVISION

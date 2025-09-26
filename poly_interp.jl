#=
Intro:
Working over polynomials F[x] over a field.
This is a tool for speeding up algebraic algorithms significantly.

Want to find a polynomial $p(x)$ given n distinct points $a_1,\ldots, a_n$ and
$y_1, \ldots, y_n$ such that $p(a_i) = y_i$. Polynomial interpolation thm says there
exists unique polynomial satisfying this with degree less than or equal to $n-1$.

We essentially get a linear system of $n$ variables with $n$ unknowns when trying to solve this.
The system is a vandermonde matrix that is invertible (ASK PAUL - vandermont?)
Vandermonde means each row is powers of some term with the first column being ones?
Determinant is called a Vandermonde determinant and has some interesting properties.
It is invertible iff all input values are distinct.

Solving using Gauss Elim takes order $n^3$ (exact formula in later lecture) time.
This is the BAD way of doing it.

Faster way requires doing order $n^2$ algorithms (with time being determined as a function
of arithmetic operations in the field).

There are two of these, Lagrange and Newton interpolation.

Lagrange Interpolation:
For Lagrange interpolation we set $L(x) := (x-a_1)\cdots(x-a_n)$ and $L_i(x) = L(x) / (x - a_i)$.
$L_i$ are called the Lagrange basis.
Then write $f := \sum_i c_i \cdot L_i$. Note, $L$ has degree $n$ and $L_i, f$ have degree $n-1$.
This gives us $f(a_i) = c_i \cdot L_i(a_i)$ and thus we want $c_i = y_i \cdot L_i(a_i)^{-1}$ 
(reciprocal in the field). Also note it's obvious we can evaluate at $a_i$ since $L_i$ doesn't have
$a_i$ as a zero (polynomials can only have n roots for integral domains).

Newton Interpolation:
For Newton interpolation we write
$f(x) = c_0 + c_1(x-a_1) + c_2(x-a_i)(x-a_2) + \ldots + c_{n-1}(x-a_1)\ldots(x-a_{n-1})$ (note - this is fairly 
similar to mixed radix form for CRT - ASK PAUL - is there a connection between these? - seems like a consequence
of CRT for polynomial vs integer rings). 
To satisfy interpolation we must have $c_0 = y_1$, $c_1 = (y_2 - c_0)/(a_2 - a_1)$ and so on.
We can compute this with Horner's form again.

Newton's interpolation is likely faster for the same reason that mixed radix form (Garner's algorithm) is.
We do want with both methods to typically convert them back to standard basis for polynomials.
Note - you can interpolate in one variable in a bivariate polynomial.

Note, there is something called the barycentric form for Lagrange's method which leads to more numerically
stable algorithms - ASK PAUL - is this worth looking into? Apparently Lagrange is useful for theoretical results,
but Newton is better in applied contexts.

Proof of uniqueness involves taking $f, g$ satisfying conditions and considering $f-g$. Then use that 
$x - a_i$ divides $f-g$ and since $x-a_i$ are coprime their product divides $f-g$ but the product has degree $n$
and $f-g$ has degree $n-1$ so $f-g = 0$.

=#

#=  Applications

Interpolation and CRT (and rational number reconstruction) are tools for speeding up algorithms - although it seems
like CRT and interpolation are two sides of the same coin.

How do we multiply two polynomials $A, B$ over integers mod $p$ with large $p$?
Instead of doing usual multiplication that we've already done (e.g. maybe some version of Karatsuba?),
we can instead think of them as polynomials and set $c(x) = a(x) b(x)$. We need to then interpolate by 
the degree of $A$ plus degree of $B$ points plus one. This avoids the need to do polynomial multiplication. Not sure
about differences in cost.

We again use Horner's method to evaluate the polynomial. This is faster for polynomials in 2 variables.

Looks like the Toom-Cook algorithms for multiplication are just to treat integers as polynomials evaluated at base
B, and then split up into $k$ smaller polynomials of about the same length, and use interpolation to multiply them
together.

For over 2 variables $x, y$, we can evaluate over $2dy + 1$ points where $dy$ is the $y$ degree of $A$, $B$. Then we
can do regular classical multiplication over the polynomials that remain to get back the product. The multiplication
classically (not using Karatsuba etc) is quadratic in $dx$ (degree of $x$). Note, we have to do this multiplication
$2dy + 1$ times since we evaluated in $y$ that many times. Then we get $2dy + 1$ points (polynomials only in $x$)
where the product evaluates to and we can reconstruct with interpolation. This takes a factor of $dy$ off the overall
computation time.

Steps:
Evaluate in $y$, $2dy + 1$ times, this is $dx+1$ polynomials (coefficients of $A$ in $x$) in both $A$ 
and $B$. Each evaluation in Horner takes $dy$ multiplications. Overall this is $2(dx+1) \cdot dy \cdot (2dy+1)$
multiplications.

Multiply each point which takes $(dx+1)^2$ multiplications, and do this $2dy+1$ times

Interpolate for $y$. The output considered as a polynomial in $x$ has degree at most $2dx$, so there are 
$2dx + 1$ polynomial coefficients in $y$. The cost of an interpolation is the degree of the polynomial squared, and 
these are in $y$ so $(2dy + 1)^2$

Overall this is order $dxdy^2 + dx^2dy$ although the coefficient of the $dy^2$ term is significantly higher.
This is a cubic algorithm - so one order of magnitude faster than the naive method. The more variables the more
speedup we get.
=#

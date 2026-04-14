#=
A lot of the intial stuff overlaps with textbook.

Total Degree
Monomial Ordering
Lexical, Graded Lexical, Grevlex ordering
ASK PAUL - Is there an ordering not satisfying well ordering that satisfies the others 
(what are we trying to exclude?)

Note - there are infinite number of monomial orderings if n >= 2. This is what makes multivariate division harder.

Assuming F is a field for now.
ALGORITHM: Division in F[x] of one polynomial by one other polynomial
Input: a, b in F[x], outputs q,r in F[x] s.t. a = bq + r

Input a, b
r <- 0
q <- 0
while a != 0 do 
    if LT(b) | LT(a) then 
        t <- LT(a) / LT(b)
        q <- q + t
        a <- a - tb
    else
        r <- r + LT(a)
        a <- a - LT(a)

output q, r

Examples in lec notes for dividing two polynomials using 2 different monomial ordering.
=#

#= Maple Polynomial Data Structures
Sum of prods rep.
Array of words - first word is header word
Header word - type of object and number of words (ie sum, prod etc)
Then we have pointer to prod arrays, followed by their coefficient
The constant is padded by a one
The product array container a header word, then pointers to variables (or e.g. functions)
followed by a word for the power.

A variable array has a name header word containing 3 things, the first is a value (if one is stored) as a pointer,
the third word is a string of the variable name (and this can be longer than 1 word if necessary).

Space: Given polynomial with t terms and n variables, memory required is <= t(2n+1+2) words
Not including memory for variables.

Multiplying monomials in this structure is very slow.
=#

#= Singular's (Oscar) Representation
Linked list representation
This package was used primarily for Groebner basis computations.
Uses a header 3 words (type, pointer to ring, pointer to polynomial first term)
Each term has a poitner to the next term (linked list), as well as the coefficient and powers
of each variable in order.

This is fast for multiplying and adding terms.
=#

#= Pari's Recursive Dense Representation
Treat it as a polynomial in the ring Z[z][y][x] (or any other ring)
so first array contains polynomial in x, with pointers to polynomials in y as the coefficients
which in turn have pointers to polynomials in z. Note - this is a dense representation.
The previous representations were sparse - both have advantages and disadvantages.
It's not truly dense since the zeros at the highest variable cut out most of the recursive data.

Here monomial multiplication is fast - just add the exponents in the tree structure.
=#

#= Maple's packed monomial arrays
We have a header word for the poly type, and then a pointer to a sequence of the variables in order.
Then we have a word for each monomial (each followed by it's coefficient). The single word is encoded
into the monomial word - where in base 10 we write e.g. 5131 where 5 is total degree then we have each degree.
They are stored using part of the word for each part of the encoding.

=#

#= Poly DAG by Monogan & Pearce
See lec notes for advantage of this for sorting and monomial comparisons.
And for additions/mults

Uses sum of prod reps if it overflows with packed being the main choice
This packing naturally works with grlex.
You could get away with not including the z index (for poly over 3 variables).

For linear polynomials maple uses sum of products representation which is better for those.
=#

#= Computing with Sparse Multivariate Polynomials
Sparse Polynomial Addition using Merging:
In Z[x] consdier h=f+g. Just do the usual merge sort.

How do multiply multivariate polynomials?
Multiply each term, then add using a merge and repeat.
Ie, f*g = (...((f1*g + f2*g) + f3*g) ... ) + fn*g = sum_{i=1,..,n} fi*g
The order of terms is maintained when multiplying by a term due to the requirement of monomial ordering.
Addition is also done by using a merge sort type of algorithm.

How do we do multivariate multiplication quickly?
NOTE - We can do Kronecker substitution then use FFT which could be quite fast!

Complexity of multiplication using merging:
Let R be an integral domain, f, g non-zero poly in R[x1,..,xn]

Let f = a1X1 + a2X2 + ... + a#X# (where Xi is a monomial and # is the number of terms in f).
Let g = b1Y1 + ... + b'Y'

and terms in f/g are sorted in a monomial ordering.

See lec notes for this stuff.....
An important optimisation within this is to turn it into a recursive algorithm by doing
getting the result of the first half of f times g plus the first half of g times f (terms of f is greater than 
terms of g) and then putting back together.

Division also works similarly. Note, speeding up division is all about understanding we are actually doing the
computation f - q*g which is mainly just a multiplication. However, we don't know q so we can't optimise in the
same way. 
We acn make it faster using:
1) Yan's geobuckets (1997) -> implemented in Singular
2) Johnson's heaps (1974) -> in Altran & Maple

Next lecture looks at 2)

We can optimise even faster via:
1) Monagon & Pearce (2008) -> in Maple 

where we can optimise via understanding that we can flip the division algorithm in some way (kind of crazy, wtf).
ASK PAUL - does he know how this one works? Where could I try to understand that?
=#

#= Using Heaps

=#


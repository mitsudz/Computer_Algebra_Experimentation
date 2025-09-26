#=
Intro:
The usual CRT approach of just modding the number by each divisor is 
actually relatively slow since modding requires divisions!

Instead, we use mixed radix form which guarantees an output which 
somehow requires less divisions? ASK PAUL - Not quite sure why this
is so much faster.
 
We can use something called Horner multiplication. -> This is just rewriting
a polynomial in a nested form then evaluating (takes out half the 
multiplications).
=#


#=
CRT for integers:

Start with $m_1,...,m_n$ and $u_1,...,u_n$. CRT says there exists unique
$u$ s.t. $u = u_i \mod m_i$ for each $i$ and $u < \prod_i m_i$

The naive way is to say that we can construct
$U = \sum_{i=1}^n u_i \prod_{j \neq i} m_j$ and $u = U mod \prod_i m_i$
This will satisfy the modulus conditions, but we have to mod by the product of divisors 
to get the second condition. We also have to work with very large numbers (potentially much 
larger than the solution).

Mixed Radix Representation: (Garner's Algorithm)
Write $u = v_1 + v2 m_1 + \ldots v_n m_1 \cdots m_{n-1}$ with each $v_i < m_i$.
Then let $s_i$ be the representation of $u$ up to the $i$'th term, and $p_i$ be the product of
$m$'s up to the $i$'th term.
Then $v_i = (u_i - s_{i-1}) * p_{i-1}^{-1} \mod m_i$ and we can compute these in a loop
without having anywhere near as many large numbers.

We use typically $m_i$ to be less than a word size and to be prime. See the lec notes about time complexity.
=#

#= 
CRT for polynomials:
We can also use this for polynomials (over integers) -> this just involves doing 
chinese remainder theorem over the coefficient for each degree. 
This also works for integer vectors componentwise.
In Maple this can be done using chrem.
We can also make the moduli polynomials for a slightly more complex set up.
Also extra details for dealing with negative moduli by using the symmetric range.
=#

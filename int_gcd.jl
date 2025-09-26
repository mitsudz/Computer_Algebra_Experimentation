# Gcd algorithms

#=
Euclid's algorithm is relatively slow since it often results in the quotient being small
for similarly sized numbers. Consecutive Fibonnaci numbers will maximise this time.

Given two random integers (what exactly do we mean by random?) the probability of getting
coprime integers is 61%. WTF - ASK PAUL - this is the reciprocal of the Riemann Zeta function
evaluated at 2???

We can do better than Euclid's algorithm. Binary GCD alg by Joseph Stein:
Use a base B of a power of 2. It is easy to test if an integer is divisible by 2 in this
base B representation. It's also relatively easy to divide by a power of 2 (O(n) work to
bit shift sufficiently - although you can also maintain a shift counter to do O(1) work ASK PAUL).

Use the following facts:
gcd(0, n) = n
gcd(a,b) = gcd(b,a)
gcd(2m, 2n) = 2 * gcd(m, n)
gcd(2m, 2n+1) = gcd(m, 2n+1)
gcd(2m+1,2n+1) = gcd(2m+1 - (2n+1), 2n+1) = gcd(m-n, 2n+1) (ensuring m >= n)

This allows us to approximately halve the size of at least one of the inputs in each case.
If one is odd and the other is even we only halve one input. If both are odd we leave the smaller
input, and more than halve the other input.

So initial guess on the time complexity of gcd(m, n) via this reduction is 
$log2 m + log2 n = log2 mn$ - was actually this squared since bitshift is linear in no. of bits

Shonhager-Strassen developed a O(m(n)logn) - this is what Maple/GMP use where 
m(n) is the multipication cost of two integers of length n.
ASK PAUL - how does this algorithm work with the matrix multiplication?
=#

""" Stein Binary GCD """
function my_gcd(m::BigInt, n::BigInt)::BigInt
    m == 0 && return n << pow
    n == 0 && return m << pow

    pow = 0
    while m > 1 && n > 1
        if iseven(m)
            m ÷= 2
            if iseven(n)
                n ÷= 2
                pow += 1
            end
            continue
        end

        # m is odd
        if iseven(n)
            n ÷= 2
            continue
        end

        # both are odd
        if m < n
            m, n = n, m
        end
        m = ((m-1)÷2) - ((n-1)÷2)
    end

    (m == 1 || n == 1) && return big(1) << pow
    return max(n, m) << pow
    # m == 0 && return n << pow
    # return m << pow
end

function my_gcd(m::Integer, n::Integer)::Integer
    return my_gcd(big(m), big(n))
end

x1 = 2596412641262175932163214713312213125964126412621759321632147133122131216321471331221312596412641262175
x2 = 21759217593216321471321632147123463324325964126216321471331221312596412641262175412621759321632147133122131
println("Number of bits is: $(Int(ceil(log2(x1+1)))), $(Int(ceil(log2(x2+1))))")
# println("My version: $(my_gcd(x1, x2))")
# println("In-built: $(gcd(x1, x2))")

# Why is this sooooo slow compared to the actual gcd implementation?
# Ask Paul for ideas, surely it's not just the fact that it's done in a higher lvl language?
# TODO - Do some sort of benchmarking and recall how to do that in Julia
#=
My function was around 45 microseconds with 33kibibytes allocated
In-built gcd was around 240 nanoseconds with only 88 bytes allocated
=#

# using Pkg
# Pkg.add("BenchmarkTools")
# Pkg.add("ProfileView")

# using BenchmarkTools
# using Profile
# using ProfileView

# show(@btime my_gcd(x1, x2))
# show(@btime gcd(x1, x2))
# @profile my_gcd(x1, x2)
# ProfileView.view()

# NOTE - Can look at @code_llvm my_gcd(x1, x2) to see the llvm code generated 
# (and @code_native for assembly code)


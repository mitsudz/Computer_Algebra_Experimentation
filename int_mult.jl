# Multiplying 2 integers of arbitrary size:
# TODO - Use maximum word size to optimise
import Base.abs
# ASK PAUL - if working on systems in Julia, will we use any languages that have a closer
# interaction with memory? Everything seems a bit too high level to optimise properly.

"""
Creates an array based representation of a (non-negative) integer `a` in a given `base``.
The smallest digit is given in position 0.
TODO - Extend this with a sign bit (or just as another type)
"""
function array_rep(a::Integer, base::Int)::Array{Int}
    if a == 0 # Edge case
        return [0]
    end
    
    @assert a > 0
    a = BigInt(a)
    n = Int(ceil(log(base, a+1))) # Number of digits needed for a
    rep = Array{Int, 1}(undef, n)
    for i = 1:n
        rep[i] = a % base
        a ÷= base
    end
    return rep
end

""" Non-Negative Integer of base B """
struct UBInt{B} <: Integer
    rep::Array{Int}

    function UBInt{B}(val::Integer) where B
        @assert B >= 2 "No unary rep"
        @assert val >= 0 "Signed not yet supported"
        new(array_rep(val, B))
    end

    function UBInt{B}(rep::Array{Int}) where B
        @assert B >= 2 "No unary rep"
        for digit in rep
            @assert digit >= 0 "Signed not supported"
            @assert digit < B "Digits must be less than $B"
        end
        new(rep)
    end
end

function Base.show(io::IO, x::UBInt{B}) where B
    for i in 1:length(x.rep)-1
        # print(Char(x.rep[end-(i-1)] - Int('0')))
        print("$(x.rep[end-(i-1)]) ")
    end
    print("$(x.rep[begin])")
end

# TODO - this wan't working, figure it out later
# @enum Sign begin # 0 can be positive or negative
#     POSITIVE
#     NEGATIVE
# end
POSITIVE = true
NEGATIVE = false

function flip(sign::Bool)::Sign
    return !sign
end

""" Integer of base B """
struct SBInt{B} <: Integer
    rep::Array{Int}
    sign::Bool

    function SBInt{B}(val::Integer) where B
        @assert B >= 2 "No unary rep"
        sign = val < 0 ? NEGATIVE : POSITIVE
        new(array_rep(abs(val), B), sign)
    end

    function SBInt{B}(rep::Array{Int}, sign::Bool) where B
        @assert B >= 2 "No unary rep"
        for digit in rep
            @assert digit >= 0 "Digits of signed integer must be non-negative"
            @assert digit < B "Digits must be less than $B"
        end
        new(rep, sign)
    end

    function SBInt{B}(a::UBInt{B}) where B
        new(a.rep, POSITIVE)
    end

    function SBInt{B}(a::UBInt{B}, sign::Bool) where B
        new(a.rep, sign)
    end
end

""" Return unsigned integer of base B """
function abs(a::SBInt{B})::UBInt{B} where B
    return UBInt{B}(a.rep)
end

function Base.show(io::IO, x::SBInt{B}) where B
    x.sign == NEGATIVE && print("-")
    show(abs(x))
end


""" Primary school multiplication of two non-negative integers represented in base B """
function primary_school_mult(a::UBInt{B}, b::UBInt{B})::UBInt{B} where B    
    res = zeros(Int, length(a.rep) + length(b.rep))
    temp = 0 # temp <= B^2 - 1 (= (B-1)^2 + 2(B-1)) - assuming for now whatever default type is large enough
    carry = 0 # carry < B
    for (i, ai) in enumerate(a.rep)
        carry = 0
        for (j, bj) in enumerate(b.rep)
            temp = res[i+j-1] + (ai*bj) + carry
            res[i+j-1], carry = temp % B, temp ÷ B
        end
        res[i+length(b.rep)] = carry
    end

    carry == 0 && pop!(res)
    return UBInt{B}(res)
end

""" Primary school multiplication of two integers represented in base B """
function primary_school_mult(a::SBInt{B}, b::SBInt{B})::SBInt{B} where B
    res = primary_school_mult(abs(a), abs(b))
    sign = b.sign == POSITIVE ? a.sign : flip(a.sign)
    return SBInt{B}(res, sign)
end

""" Karatsuba's algorithm
Assume for simplicity the length of both numbers is the same and a power of 2, $2^k$.
We write $a = a1|a2$, $b = b1|b2$ where $a1,a2,b1,b2$ have length $n/2$.
Then, $ab = (a1 \cdot B^{2_{k-1}} + a2)(b1 \cdot B^{2_{k-1}} + b2)$
$ = a1b1 \cdot B^{2_{k}} + (a1b2 + a2b1) \cdot B^{2_{k-1}} + a2b2$
This is 4 multiplications of integers half as long, plus multiplying by powers of B (which is just shifting).
We do this recursively until one number is 0, or they are 1 digit long.

We can make this faster by using the following trick:
$a1b2 + a2b1 = (a1-a2)(b2-b1) + a1b1 + a2b2$
This has 3 multiplications, but we are already doing the last two and addition is linear 
and is thus faster. The coefficient of the linear work for addition/subtraction increases,
but the recursive splitting at each level decreases to 3 and so we get an upper bound on time
complexity of $n^{\log_2(3)}$ rather than $n^2$ as in the above.

Overall, we do:
$ab = a1b1 \cdot B^{2_{k}} + ((a1-a2)(b2-b1) + a1b1 + a2b2) \cdot B^{2_{k-1}} + a2b2$

We use B = word size (ie size of Int)

There are numerous better complexity algorithms including Toome's various ones,
Strassen and S... and Harvey... which get it down to O(nlogn)
"""
function __karatsuba_rec(a::Array{Int}, b::Array{Int}, base::Integer)::Array{Int}
    if all(iszero, a) || all(iszero, b)
        return UBInt(0)
    end

    @assert false "base case for length 1 arrays"

    halfN = length(a.rep) ÷ 2
    a1 = a[1+halfN : end]
    a2 = a[begin : halfN]
    b1 = b[1+halfN : end]
    b2 = b[begin : halfN]

    a1b1 = __karatsuba_rec(a1, b1, base)
    a2b2 = __karatsuba_rec(a2, b2, base)
    cross_term = 0
    @assert false "Not yet implemented"
end

function karatsuba(a::UBInt{B}, b::UBInt{B})::UBInt{B} where B
    # Pad to the next power of 2 for simplicity
    n = max(nextpow(2, length(a.rep)), nextpow(2, length(b.rep))
    if length(a.rep) < n
        m = length(a.rep)
        resize!(a.rep, n)
        fill!(view(a.rep, m+1: n), 0) # Should use zero(Int)
    end
    if length(b.rep) < n
        m = length(b.rep)
        resize!(b.rep, n)
        fill!(view(b.rep, m+1: n), 0) # Should use zero(Int)
    end

    return UBInt(__karatsuba_rec(a.rep, b.rep, B), B)
end

# Testing
base = 10

x1 = UBInt{base}(15)
x2 = UBInt{base}(12)
println(x1, x2)
println(primary_school_mult(x1, x2))

x1 = SBInt{base}(-6)
x2 = SBInt{base}(11)
println(x1, x2)
println(primary_school_mult(x1, x2))

x1 = UBInt{base}(15)
x2 = UBInt{base}(12)
println(x1, x2)
println(karatsuba(x1, x2))

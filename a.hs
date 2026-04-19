Invariant:
normTR x0 {x1,...,xn} => x0^2 + ||(x1,...,xn)||^2
                       = x0^2 + x1^2 + ... + xn^2

normTR ans (x1:xs) = normTR sqrt(ans^2 + x1^2) xs

LHS Invariant: ans^2 + x1^2 + x2^2 + ... + xn^2
RHS Invariant: (sqrt(ans^2 + x1^2))^2 + x2^2 + x3^2 + ... + xn^2
             = (ans^2 + x1^2) + x2^2 + ... + xn^2
             = ans^2 + x1^2 + x2^2 + ... +xn^2
             = LHS Invariant

normTR ans [] = ans
LHS Invariant: ans^2
RHS = sqrt(LHS Invariant)

normTR 567 [2, 41, 1, 5, 6] = ... = normTR ans []
ans^2 = 567^2 + ... + 6^2

norm xs = normTR 0 xs = ... = normTR ans [] = ans
=> 0 + x1^2 + x2^2 + ... + xn^2 = ans^2
=> sqrt(0 + x1^2 + x2^2 + ... + xn^2) = sqrt(ans^2) = ans

|| xs || = sqrt(x1^2 + x2^2 + ... + xn^2)

normTR 2 [1, 2, 5, 9] => 2^2 + 1^2 + 2^2 + 5^2 + 9^2 = 115
= normTR sqrt(5) [2, 5, 9] => 5 + 2^2 + 5^2 + 9^2 = 115

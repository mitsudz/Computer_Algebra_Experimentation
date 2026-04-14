#=
Use a heap to maintain the next largest terms in the product that 
potentially have the same degree. Do this by ordering all the products 
(implicitly) and always adding either the next element and sometimes
the start of the next row.

At any point in the algorithm you only ever have at most one element
from every row of the polynomial product.

Cost becomes O(l*m*log(m)) where m, l are the number of non-zero terms
resp. Should choose m to be the smaller of the two. Naive multiplication
is O(m^2*l).
=#

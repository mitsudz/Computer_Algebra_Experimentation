#= How to optimise in Julia
Inbounds
Multiprocessing/threading
Hoisting - moving fixed memory access outside of loops
Value Types - give info to compiler to do optimisation whilst allowing more generic code - don't need to name this!
Note - by giving value types for fixed length arrays this does optimise by doing loop unrolling without us having to 
manually do so.
SVector macro for static arrays
Moving up base case of recursions
LoopVectorisation.jl has @Turbo (note we also have a @SIMD macro somewhere)! @turbo is a bit fudged atm - doesn't work
perfectly to vectorise the code. But when it does work it is hella fast.

=#

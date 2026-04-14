#############################################################################
#############################################################################
#
# This file implements the Heap data structure 
#                                                                               
#############################################################################
#############################################################################

import Base: length, push!, pop!, show, peek, isempty


""" 
A struct to represent a basic array-based max-heap.

T -> can be any type which implements comparison (e.g., integers, reals etc.)

The underlying array representation can be accessed via h.data where h::Heap.

HEAP ORDER:
    See the linked article in the instructions for Task 3 to understand heap order.
    In particular, we must have that given a vector `vec`, `vec` is in (max) heap order if
    for any index `i`, we have `vec[i] ≥ vec[2 * i]` and `vec[i] ≥ vec[2 * i + 1]` (assuming
    `length(vec) ≥ 2 * i` and `length(vec) ≥ 2 * i + 1`).
    
    Note: for terms in sparse polynomials stored in a heap, the heap order inequality should be strict.
"""
struct Heap{T}
    data::Vector{T}

    """ Empty heap constructor """
    Heap{T}() where {T} = new(T[])

    """ 
    Inner constructor to create a heap from an array.

    WARNING:
        Does NOT check whether the heap property is maintained and does not enforce this.
        You must ensure the heap property is satisfied (e.g., the array is sorted in descending order)
        before using this constructor.  

        Do NOT use this constructor to create a heap from an arbitrary vector - use an outer constructor.
    """
    Heap{T}(v::Vector{T}) where {T} = new(v) 
end

#####################
# Heap constructors #
#####################

""" Convenience function to create an empty heap containing elements of type T. """
Heap(::Type{T}) where {T} = Heap{T}()

"""
Given a vector, turns it into a heap (establishing heap order).
Stores a copy of the vector (the original is unchanged).

This function is safe to use for an arbitrary vector.
"""
function Heap(v::Vector{T}) where {T}
    return Heap!(deepcopy(v))
end

"""
Given a vector, turns it into a heap (establishing heap order).
This function stores the original vector (i.e., the original vector may be modified through 
use/construction of the heap).

This function is safe to use for an arbitrary vector.
"""
function Heap!(v::Vector{T}) where {T}
    h = Heap{T}(v)
    _heapify!(h)
    return h
end

"""
Maps elements (non-destructively) of the heap via an order preserving function.
The original heap is left unchanged.

WARNING: 
    This will break the heap if the `order_pres_func` is not order preserving.
    An order preserving function `f` will maintain the order of any two inputs.
    I.e., if `a > b` then `f(a) > f(b)`
"""
function map_heap(h::Heap{T}, order_pres_func::Function) where {T}
    Heap{T}(order_pres_func.(h.data))
end

"""
Maps elements (destructively) of the heap via an order preserving function.
The original heap is modified in-place.

WARNING: See warning above for map_heap
"""
function map_heap!(h::Heap{T}, order_pres_func::Function) where {T}
    h.data .= order_pres_func.(h.data)
    return h
end

###########
# Display #
###########

show(io::IO, h::Heap) = print(io, "Heap(", h.data, ")")

#####################
# Querying the heap #
#####################

""" Get the number of elements in the heap. """
length(h::Heap) = length(h.data)

""" Check if the heap is empty """
isempty(h::Heap)::Bool = length(h) == 0

""" 
Push a new element onto the heap.
The heap will automatically handle maintaining heap order
"""
function push!(h::Heap{T}, x::T) where {T}
    push!(h.data, x)
    _up_heap!(h, length(h.data))
end

""" 
Pop (remove and return) the maximum element from the heap.
The heap will automatically handle maintaining heap order.

Note: This function will error if called when the heap is empty - check first with isempty(h)
"""
function pop!(h::Heap)
    isempty(h) && error("Heap is empty, cannot pop!")
    length(h) == 1 && return pop!(h.data)

    # Swap first/last elements
    max_val = h.data[1]
    h.data[1] = pop!(h.data)
    
    # Restore heap order
    if !isempty(h)
        _down_heap!(h, 1)
    end
    
    return max_val
end

"""
Pops all elements from the heap, and returns a vector (in sorted order).

Note: the original heap will be empty at the end of this operation.
"""
function popall!(h::Heap{T})::Vector{T} where {T}
    n = length(h)
    v = Vector{T}()
    for _ in 1:n
        push!(v, pop!(h))
    end
    @assert isempty(h)
    return v
end

"""
Pops all elements from the heap, and returns a vector (in sorted order).

Note: the original heap will unaffected.
"""
function popall(h::Heap{T})::Vector{T} where {T}
    return popall!(deepcopy(h))
end

"""
Returns the maximum element from the heap without removing it.

Note: This function will error if called when the heap is empty - check first with isempty(h)
"""
function peek(h::Heap) 
    isempty(h) && error("Heap is empty, cannot peek")
    return h.data[1]
end

####################################################
# Internal helper functions to maintain heap order #
####################################################

"""
Reorders a vector in place to satisfy the heap property.
"""
function _heapify!(h::Heap)
    n = length(h.data)
    for i in div(n, 2):-1:1
        _down_heap!(h, i)
    end
end

""" Internal helper function to restore heap order upwards """
function _up_heap!(h::Heap, i::Integer)
    parent_index = div(i, 2)
    while i > 1 && h.data[i] > h.data[parent_index]
        h.data[i], h.data[parent_index] = h.data[parent_index], h.data[i]
        i = parent_index
        parent_index = div(i, 2)
    end
end

""" Internal helper function to restore heap order downwards """
function _down_heap!(h::Heap, i::Integer)
    n = length(h.data)
    while (left_child = 2 * i) <= n
        j = left_child
        if (right_child = 2 * i + 1) <= n && h.data[right_child] > h.data[left_child]
            j = right_child
        end

        h.data[i] > h.data[j] && break
        h.data[i], h.data[j] = h.data[j], h.data[i]
        i = j
    end
end

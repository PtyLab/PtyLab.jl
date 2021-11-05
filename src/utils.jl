function circ(shape, xi, radius)
    rr2 = IndexFunArrays.rr2(shape, scale=ScaMid, offset=CtrMid) .* xi[end]^2
    return rr2 .< radius^2
end


"""
    randpermPositions(arr::AbstractArray{T, 2}) where T

Random shuffle the columns of a matrix.


 # Examples
```julia
julia> randpermPositions(([1 2 3 4 5; 1 2 3 4 5]))
2×5 Matrix{Int64}:
 5  4  2  1  3
 5  4  2  1  3

julia> randpermPositions(([1 2 3 4 5; 1 2 3 4 5]))
2×5 Matrix{Int64}:
 2  5  1  4  3
 2  5  1  4  3
```
"""
function randpermPositions(arr::AbstractArray{T, 2}) where T
    out = similar(arr) 
    
    # random permutation of the numbers
    rperm = randperm(size(arr, 2))
    for (i_new, i_rand) in enumerate(rperm)
        # no copy with view
        out[:, i_new] .= view(arr, :, i_rand)
    end
    return out, rperm 
end

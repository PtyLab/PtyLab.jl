export Fraunhofer

"""
    Fraunhofer(arr, rec, params; shift=false, dims=(1,2))

Returns two functions `object2detector, detector2object` which can propagate `asw` with Fraunhofer assumption
efficiently between object and detector back and forth


 ## Examples
```julia-repl
julia> 
```
"""
function Fraunhofer(arr, rec, params; shift=false, dims=(1,2))
    p = plan_fft(arr, dims)
  
    # create object to detector 
    o2d(x) = p * x
    # same with shift
    o2ds(x) = fftshift(p * ifftshift(x, dims), dims)

    # create detector to object 
    d2o(x) = inv(p) * x 
    d2os(x) = fftshift(inv(p) * ifftshift(x, dims), dims)

    if shift
        return o2ds, d2os
    else
        return o2d, d2o
    end
end

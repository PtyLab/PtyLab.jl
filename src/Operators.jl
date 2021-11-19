export Fraunhofer

"""
    Fraunhofer(arr; fftshift=false, dims=(1,2), FFTW_flags=FFTW.MEASURE)

Returns two functions `object2detector, detector2object` which can propagate `asw` with Fraunhofer assumption
efficiently between object and detector back and forth
Currently uses `plan_fft!` for in-place ffts. `FFTW_flags` is only passed to `plan_fft` if the array is a CPU array.

 ## Examples
 `arr` should be something like the ESW `rec.object[1:rec.Np, 1:rec.Np, ..] .* rec.probe`.
```julia-repl
julia> o2d, d2o = Fraunhofer(arr)
```
"""
function Fraunhofer(arr::T; fftshift=false, dims=(1,2), FFTW_flags=FFTW.PATIENT) where T

    # planning can overwritte array sometimes!
    # only FFTW allows different flags
    p = if T <: CuArray
            plan_fft!(similar(arr), dims)
        else
            plan_fft!(similar(arr), dims, flags=FFTW_flags)
    end

    # those let blocks are required because of 
    # https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured


    # create object to detector 
    o2d! = let p=p
        o2d!(x) = p * x
    end
    # same with shift
    o2ds = let p=p
        o2ds(x) = fftshift(p * x, dims)
    end

    # create detector to object without shifts 
    d2o! = let p_inv=inv(p)
        d2o!(x) = p_inv * x
    end

    d2os = let p_inv=inv(p)
        d2os(x) = p_inv * ifftshift(x, dims)
    end


    # return correct functions
    if fftshift 
        return o2ds, d2os
    else
        return o2d!, d2o!
    end
end

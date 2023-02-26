export Fraunhofer, ASPW


"""
    _plan_fft!_cuda_FFTW(T, arr, dims, FFTW_flags)
    
"""
function _plan_fft!_cuda_FFTW(T, arr, dims, FFTW_flags)
    # planning can overwritte array sometimes!
    # only FFTW allows different flags
    p = if T <: CuArray
            plan_fft!(similar(arr), dims)
        else
            plan_fft!(similar(arr), dims, flags=FFTW_flags)
    end
    return p
end



"""
    Fraunhofer(arr; fftshiftSwitch=false, dims=(1,2), FFTW_flags=FFTW.MEASURE)

Returns two functions `object2detector, detector2object` which can propagate `asw` with Fraunhofer assumption
efficiently between object and detector back and forth
Currently uses `plan_fft!` for in-place ffts. `FFTW_flags` is only passed to `plan_fft` if the array is a CPU array.

 ## Examples
 `arr` should be something like the ESW `rec.object[1:rec.Np, 1:rec.Np, ..] .* rec.probe`.
```julia-repl
julia> object2detector, detector2object = Fraunhofer(arr)
```
"""
function Fraunhofer(arr::T; fftshiftSwitch=false, dims=(1,2), FFTW_flags=FFTW.MEASURE) where T


    ss = sqrt(size(arr, 1) * size(arr, 2))
    # those let blocks are required because of 
    # https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured
    p = _plan_fft!_cuda_FFTW(T, arr, dims, FFTW_flags)

    # create object to detector 
    object2detector! = let p=p
        function object2detector!(x)
            px = p * x
            px ./= ss
        end
    end
    # same with shift
    object2detectors = let p=p
        function object2detectors(x)
            px = fftshift(p * x, dims)
            px ./= ss
        end
    end

    # create detector to object without shifts 
    detector2object! = let p_inv=inv(p)
        function detector2object!(x)
            px = p_inv * x
            px .*= ss 
        end
    end

    detector2objects = let p_inv=inv(p)
        function detector2objects(x)
            px = p_inv * ifftshift(x, dims)
            px .*= ss 
        end
    end


    # return correct functions
    if fftshiftSwitch 
        return object2detectors, detector2objects
    else
        return object2detector!, detector2object!
    end
end

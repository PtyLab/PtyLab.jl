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
    Fraunhofer(arr; fftshiftFlag=false, dims=(1,2), FFTW_flags=FFTW.MEASURE)

Returns two functions `object2detector, detector2object` which can propagate `asw` with Fraunhofer assumption
efficiently between object and detector back and forth
Currently uses `plan_fft!` for in-place ffts. `FFTW_flags` is only passed to `plan_fft` if the array is a CPU array.

 ## Examples
 `arr` should be something like the ESW `rec.object[1:rec.Np, 1:rec.Np, ..] .* rec.probe`.
```julia-repl
julia> o2d, d2o = Fraunhofer(arr)
```
"""
function Fraunhofer(arr::T; fftshiftFlag=false, dims=(1,2), FFTW_flags=FFTW.MEASURE) where T


    # those let blocks are required because of 
    # https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured
    p = _plan_fft!_cuda_FFTW(T, arr, dims, FFTW_flags)

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
    if fftshiftFlag 
        return o2ds, d2os
    else
        return o2d!, d2o!
    end
end

"""
    ASPW(u::A, z::T, wavelength::T, L::T; dims=(1,2), FFTW_flags=FFTW.MEASURE)


Create a method `d2o(u, H=H)` which can efficienty propagate a field with
the angular spectrum method.

# TODO
* untested
"""
function ASPW(u::A, z::T, wavelength::T, L::T; dims=(1,2), FFTW_flags=FFTW.MEASURE) where {A, T}
    p = _plan_fft!_cuda_FFTW(A, u, dims, FFTW_flags)

    d2o = let 
        H = ASPW_kernel(u, z, wavelength, L) 
        function d2o(u, H=H)
            # fourier transforms are in-place
            U = p * u
            u = inv(p) * (U .* H)
            return u 
        end
    end
    return d2o
end

"""
    ASPW_kernel(u::A, z::T, wavelength::T, L::T) where {A, T}

Calculates the kernel which is multiplied in Fourier space for the ASPW propagation.
It is centered around `(1,1)`.
"""
function ASPW_kernel(u::A, z::T, wavelength::T, L::T) where {A, T}
    k = T(2π) / wavelength
    N = size(u, 1)

    Fx = similar(u, real(eltype(u)), N)
    Fx .= fftfreq(N, N / L)
    Fy = Fx'
    f_max = L / (wavelength * sqrt(L^2 + 4 * z^2))

    W = ifftshift(circ.(Fx, Fy, 2 * f_max))
    H = W .* cis.(k * z * sqrt.(0im .+ 1 .- (Fx .* wavelength).^2 .- (Fy .* wavelength).^2)) 
    
    return H
end

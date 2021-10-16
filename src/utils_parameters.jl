import Base.@kwdef

export init_Parameters

"""
    Parameters{T}

Struct storing all parameters.
"""
@kwdef struct Parameters{T}
    wavelength::T
    # probe sampling
    dxp::T
    Np::Int
    xp::Vector{T}
    Lp::T
    zp::Union{T, Nothing}
    # object sampling
    dxo::T
    No::Int
    xo::Vector{T}
    Lo::T
    zo::T
    # detector sampling
    dxd::T
    Nd::Int
    xd::Vector{T}
    Ld::T
    # the general data type which is enforced
    dtype::Type{T}
end

"""
    init_Parameters(<kwargs>)

Function to return parameter object storing
all meta information needed for reconstruction.
Note that can access all those parameters but many
of them are connected hence to fill the struct we only need a few of
them.

# Parameters
 
## Physical Properties
* `wavelength` of the laser

## General Properties

* `dtype`: datatype used for reconstructions. `Float32` is usually faster, especially on GPUs.

## Probe
* `dxp`: pixel size
* `Np`: number of pixels
* `xp`: 1D coordinates 
* `Lp`: field of view probe (field size)
* `zp`: distance to next plane

## Object 
* `dxo`: pixel size
* `No`: number of pixels
* `xo`: 1D coordinates 
* `Lo`: field of view probe (field size)
* `zo`: distance to next plane

## Detector
* `dxd`: pixel size
* `Nd`: number of pixels
* `xd`: 1D coordinates 
* `Ld`: field of view probe (field size)

"""
function init_Parameters(;
        wavelength=DEFAULT_WAVELENGTH, 
        No=2^7,
        Nd=2^9,
        dxd=4.5e-6,
        zo=50e-3,
        dtype=Float32,
        zp=nothing
        )

    Ld = Nd * dxd
    # probe stuff
    dxp = wavelength * zo / Ld       
    Np = Nd
    Lp = dxp*Np
    xp = FourierTools.fftpos(Lp, Np) 
    # object stuff
    dxo = dxp
    Lo = dxo * No
    xo = FourierTools.fftpos(Lo, No) 
    xd = FourierTools.fftpos(Ld, Nd) 

    _dtype(x) = isnothing(x) ? x : dtype(x)
    
    Parameters(
        wavelength=_dtype(wavelength),
        # probe sampling
        dxp=_dtype(dxp),
        Np=Np,
        xp=_dtype.(xp),
        Lp=_dtype(Lp),
        zp=_dtype(zp),
        # object sampling,
        dxo=_dtype(dxo),
        No=No,
        xo=_dtype.(xo),
        Lo=_dtype(Lo),
        zo=_dtype(zo),
        # detector sampling,
        dxd=_dtype(dxd),
        Nd=Nd,
        xd=_dtype.(xd),
        Ld=_dtype(Ld),
        # the general data type which is enforced,
        dtype=dtype)
end



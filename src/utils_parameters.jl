import Base.@kwdef

export Parameters

@kwdef struct Parameters{T}
    wavelength::Union{T, Nothing}
    # probe sampling
    dxp::T
    Np::Int
    xp::Vector{T}
    Lp::T
    zp::Union{T, Nothing}
    # object sampling
    dxo::T
    No::Int
    xo::T
    Lo::T
    zo::T
    # detector sampling
    dxd::T
    Nd::Int
    xd::Vector{T}
    Ld::T
    # the general data type which is enforced
    dtype::Type{T}
    
    # function Parameters(;
    #     wavelength = 632.8e-9,
    #     # probe
    #     dxp=nothing,
    #     Np=2^7,
    #     Lp=nothing,
    #     zp=nothing,
    #     # object
    #     dxo=nothing,
    #     No=2^7,
    #     Lo=nothing,
    #     zo=nothing,
    #     # detector
    #     dxd=dxp,
    #     Nd=2^9,
    #     Ld=nothing,
    #     )

end

function Parameters(;
        wavelength=DEFAULT_WAVELENGTH, 
        No=2^7,
        Nd=2^9,
        dxd=4.5e-6,
        zo=50e-3,
        dtype=Float32,
        zp=nothing
        )

    # @show "lol"
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
    
    Parameters(;
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
        Ld=Ld,
        # the general data type which is enforced,
        dtype=dtype
    )
end



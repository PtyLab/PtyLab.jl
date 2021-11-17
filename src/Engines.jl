export reconstruct
export ePIE

abstract type Engines end

"""
    @with_kw mutable struct ePIE{T} <: Engines

`ePIE` is a struct to store important parameters for the `ePIE::Engine` 
"""
@with_kw mutable struct ePIE{T} <: Engines
    betaProbe::T=0.25f0
    betaObject::T=0.25f0
    numIterations::Int=50
end

"""
    IntensityProjection(rec::ReconstructionCPM{T}, params::Params)

Returns a function which does the `intensityProjection`.

 # Examples
```julia
julia> intensityProjection = IntensityProjection(rec, params)

julia> eswUpdate = intensityProjection(esw, Imeasured)
```
"""
function IntensityProjection(rec::ReconstructionCPM{T}, params::Params) where T
    # create an efficient propagator function
    esw_temp = rec.object[1:rec.Np, 1:rec.Np, ..] .* rec.probe
    object2detector, detector2object = params.propagatorType(esw_temp, params)
       
    
    @warn "gimmel is currently estimated as `100 * eps($T)`"
    gimmel = 100 * eps(T)
    f! = let    intensityConstraint = params.intensityConstraint
                object2detector = object2detector
                detector2object = detector2object
        function f(esw, Imeasured)
            ESW = object2detector(esw)
            Iestimated = let
                if intensityConstraint === IntensityConstraintStandard
                    # sum over the last three channels.
                    # @tullio Iestimated[i, j] := abs2(ESW[i,j,k,s1,s2,s3]) 
                    # that is currently faster than @tullio
                    @tturbo view(sum(abs2, ESW, dims=(3, 4, 5, 6)), :, :, 1,1,1,1)
                else
                    error("Unknown intensityConstraint $intensityConstraint")
                end
            end

            frac = let 
                if intensityConstraint === IntensityConstraintStandard 
                    @tullio frac[a1,a2] := sqrt(Imeasured[a1,a2] / (Iestimated[a1,a2] + gimmel))
                    # frac = sqrt.(Imeasured ./ (Iestimated .+ gimmel))
                else
                    error("Unknown intensityConstraint")
                end
            end

            # update ESW
            ESW .*= frac

            # back to detector, memory free due to plan_fft!
            eswUpdate = detector2object(ESW)
            return eswUpdate 
        end
    end

    return f! 
end

"""
    probeUpdate(engine::ePIE{T}, objectPatch, probe, DELTA) where T

Returns the `newProbe`.

TODO memory improvements possible

"""
function probeUpdate(engine::ePIE{T}, objectPatch, probe, DELTA) where T
    fracProbe = conj.(objectPatch) ./ maximum(sum(abs2.(objectPatch), dims=3:ndims(objectPatch)))
    newProbe = probe .+ engine.betaProbe .* sum(fracProbe .* DELTA, dims=(3, 5, 6))

    return newProbe
end

"""
    probeUpdate(engine::ePIE{T}, objectPatch, probe, DELTA) where T

Returns the `newobjectPatch`.
TODO memory improvements possible
"""
function objectPatchUpdate(engine::ePIE{T}, objectPatch, probe, DELTA) where T
    fracObject = conj.(probe) ./ maximum(sum(abs2.(probe), dims=3:ndims(probe)))
    newObject = objectPatch .+ engine.betaObject .* sum(fracObject .* DELTA, dims=(3, 4, 6))

    return newObject
end

"""
    reconstruct(engine::ePIE{T}, params::Params, rec::ReconstructionCPM{T}) where T 

Reconstruct a CPM dataset.
"""
function reconstruct(engine::ePIE{T}, params::Params, rec::ReconstructionCPM{T}) where T 
    # calculate the positions since rec.positions is in fact a function!
    positions = rec.positions

    # create intensityProjection function 
    intensityProjection! = IntensityProjection(rec, params)

    # get two function which maybe shift depending on the flag             \
    maybe_fftshift, maybe_ifftshift = getMaybeFftshifts(! params.fftshiftFlag) 

    # alias and shift maybe
    object = rec.object
    probe = rec.probe
    ptychogram = maybe_ifftshift(rec.ptychogram)

    Np = rec.Np

    # copy them!
    # those are buffers
    oldProbe = copy(probe)
    oldObjectPatch = object[1:Np, 1:Np, ..]
    esw = object[1:Np, 1:Np, ..] .* probe
    
    @showprogress for loop in 1:engine.numIterations
        for positionIndex in randperm(size(positions, 2))
            # row, col does not work!
            col, row = positions[:, positionIndex] 
                
            sy = row:(row + Np - 1)
            sx = col:(col + Np - 1)
            # no copy necessary -> because esw immediately calculated
            objectPatch = view(object, sy, sx, ..)

            # save old state since we need that later for probeUpdate and objectPatchUpdate
            oldObjectPatch .= objectPatch
            oldProbe .= probe

            # exit surface wave,
            # tullio seems to be slower on such simple operations
            # @tullio esw[i1,i2,i3,i4,i5,i6] = objectPatch[i1,i2,i3,i4,1,i6] .* probe[i1,i2,i3,1,i5,i6]
            esw = objectPatch .* probe

            # already store esw in DELTA, since intensityProjection is going to change esw
            DELTA = -1 .* esw
            eswUpdate = intensityProjection!(esw, view(ptychogram, :, :, positionIndex))
            
             # difference term
            DELTA .+= eswUpdate

            # update newProbe und newObjectPatch
            newProbe = probeUpdate(engine, oldObjectPatch, oldProbe, DELTA) 
            newObjectPatch = objectPatchUpdate(engine, oldObjectPatch, oldProbe, DELTA) 
            probe .= newProbe
            object[sy, sx, ..] .= newObjectPatch
        end 
    end
    return probe, object
end

export reconstruct
export ePIE

abstract type Engines end


@with_kw mutable struct ePIE{T} <: Engines
    betaProbe::T=0.25f0
    betaObject::T=0.25f0
    numIterations::Int=50
end


function IntensityProjection(rec::ReconstructionCPM{T}, params::Params) where T
    # create an efficient propagator function
    esw_temp = rec.object[1:rec.Np, 1:rec.Np, ..] .* rec.probe
    object2detector, detector2object = params.propagatorType(esw_temp, rec, params)
       
    
    @warn "gimmel is currently estimated as `100 * eps($T)`"
    gimmel = 100 * eps(T)
    function f(esw, Imeasured)
        ESW = object2detector(esw)
        Iestimated = let
            if params.intensityConstraint === IntensityConstraintStandard
                # sum over the last three channels.
                sum(abs2, ESW, dims=(4, 5, 6))
            else
                error("Unknown intensityConstraint")
            end
        end

        frac = let 
            if params.intensityConstraint === IntensityConstraintStandard 
                frac = sqrt.(Imeasured ./ (Iestimated .+ gimmel))
            else
                error("Unknown intensityConstraint")
            end
        end

        # update ESW
        ESW .*= frac

        # back to detector
        eswUpdate = detector2object(ESW)
        return eswUpdate 
    end

    return f 
end


function probeObjectPatchUpdate(engine::ePIE{T}, objectPatch, probe, DELTA) where T
    fracObject = conj.(probe) ./ maximum(sum(abs2.(probe), dims=3:ndims(probe)))
    fracProbe = conj.(objectPatch) ./ maximum(sum(abs2.(objectPatch), dims=3:ndims(objectPatch)))

    newProbe = probe .+ engine.betaProbe .* sum(fracProbe .* DELTA, dims=(3, 5, 6))
    newObject = objectPatch .+ engine.betaObject .* sum(fracObject .* DELTA, dims=(3, 4, 6))

    return newProbe, newObject
end

function reconstruct(engine::ePIE{T}, params::Params, rec::ReconstructionCPM{T}) where T 
    # calculate the positions
    positions = rec.positions


    # create a 
    intensityProjection = IntensityProjection(rec, params)
    object = rec.object
    probe = rec.probe
    ptychogram = rec.ptychogram
    @showprogress for loop in 1:engine.numIterations
        Np = rec.Np
        for positionIndex in randperm(size(positions, 2))
            # row, col does not work!
            col, row= positions[:, positionIndex] 
                
            sy = row:(row + Np - 1)
            sx = col:(col + Np - 1)
            # using EllipsisNotation (..)
            # this line already copies the data!
            objectPatch = object[sy, sx, ..]

            # exit surface wave
            esw = objectPatch .* probe

            eswUpdate = intensityProjection(esw, view(ptychogram, :, :, positionIndex))
                
             # difference term
            DELTA = eswUpdate .- esw

            # update newProbe und newObjectPatch
            newProbe, newObjectPatch = probeObjectPatchUpdate(engine, objectPatch, probe, DELTA) 
            probe .= newProbe
            object[sy, sx, ..] .= newObjectPatch 
        end 
        # break
    end
    return rec.probe, rec.object
end

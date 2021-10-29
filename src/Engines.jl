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
    object2detector, detector2object = params.propagatorType(rec.object, rec, params)
       
    
    @warn "gimmel is currently estimated as `100 * eps($T)`"
    gimmel = 100 * eps(T)
    function f(esw, Imeasured)
        ESW = object2detector(esw)
        Iestimated = let
            if params.intensityConstraint === IntensityConstraintStandard
                sum(abs2, ESW, dims=(1, 2))
            end
        end

        frac = let 
            if params.intensityConstraint === IntensityConstraintStandard 
                frac = sqrt.(Imeasured ./ (Iestimated .+ gimmel))
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

    # def objectPatchUpdate(self, objectPatch: np.ndarray, DELTA: np.ndarray):
    #     """
    #     Todo add docstring
    #     :param objectPatch:
    #     :param DELTA:
    #     :return:
    #     """
    #     # find out which array module to use, numpy or cupy (or other...)
    #     xp = getArrayModule(objectPatch)

    #     frac = self.reconstruction.probe.conj() / xp.max(xp.sum(xp.abs(self.reconstruction.probe) ** 2, axis=(0, 1, 2, 3)))
    #     return objectPatch + self.betaObject * xp.sum(frac * DELTA, axis=(0,2,3), keepdims=True)

       
    # def probeUpdate(self, objectPatch: np.ndarray, DELTA: np.ndarray):
    #     """
    #     Todo add docstring
    #     :param objectPatch:
    #     :param DELTA:
    #     :return:
    #     """
    #     # find out which array module to use, numpy or cupy (or other...)
    #     xp = getArrayModule(objectPatch)
    #     frac = objectPatch.conj() / xp.max(xp.sum(xp.abs(objectPatch) ** 2, axis=(0,1,2,3)))
    #     r = self.reconstruction.probe + self.betaProbe * xp.sum(frac * DELTA, axis=(0, 1, 3), keepdims=True)
    #     return r




function probeObjectPatchUpdate!(engine::ePIE{T}, objectPatch, probe, DELTA) where T
    fracProbe = conj.(probe) ./ maximum(sum(abs2.(probe), dims=3:ndims(probe)))
    fracObject = conj.(objectPatch) ./ maximum(sum(abs2.(objectPatch), dims=3:ndims(objectPatch)))

    newProbe .= probe .+ engine.betaProbe .* sum(frac .* DELTA, dims=(2, 3, 5))
    newObject .= objectPatch .+ engine.betaObject .* sum(frac .* DELTA, dims=(2,4,5))

    return newProbe, newObject
end

function reconstruct(engine::ePIE{T}, params::Params, rec::ReconstructionCPM{T}) where T 
    # calculate the positions
    positions = rec.positions


    # create a 
    intensityProjection = IntensityProjection(rec, params)
    for loop in 1:engine.numIterations
        Np = rec.Np
        @warn "Order is not randomized yet"
        for positionIndex = 1:rec.numFrames
            @show positions
            row, col = positions[positionIndex] 
                
            sy = row:(row + Np)
            sx = col:(col+ Np)
            # using EllipsisNotation
            # this line already copies the data!
            objectPatch = rec.object[sy, sx, ..]

            # exit surface wave
            esw = objectPatch .* rec.probe

            eswUpdate = intensityProjection(esw, rec.ptychogram[positionIndex])
                
             # difference term
            DELTA .-= esw

            probeObjectPatchUpdate!(engine, objectPatch, rec.probe, DELTA) 
        end 
    end
end

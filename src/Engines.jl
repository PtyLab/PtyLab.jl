
abstract type Engines{T} end


@struct mutable struct ePIE{T} <: Engines{T}
    betaProbe::T=T(0.25)
    betaObject::T=T(0.25)
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
    
        # frac = sqrt.(self.reconstruction.Imeasured / (self.reconstruction.Iestimated + gimmel))
    end


end


function reconstruct(engine::Engines{T}, params::Params, rec::ReconstructionCPM{T}) where T 
    # calculate the positions
    positions = rec.positions


    # create a 
    intensityProjection = IntensityProjection(rec)  
    for loop in 1:engine.numIterations
        Np = rec.Np
        @warn "Order is not randomized yet"
        for positionIndex = 1:rec.numFrames
            row, col = rec.positions[positionIndex] 
                
            sy = row:(row + Np)
            sx = col:(col+ Np)
            # using EllipsisNotation
            # this line already copies the data!
            objectPatch = rec.object[sy, sx, ..]

            # exit surface wave
            esw = objectPatch .* rec.probe

            # propagate to camera and 
            
                # make exit surface wave
                self.reconstruction.esw = objectPatch * self.reconstruction.probe

                # propagate to camera, intensityProjection, propagate back to object
                self.intensityProjection(positionIndex)

                # difference term
                DELTA = self.reconstruction.eswUpdate - self.reconstruction.esw

                # object update
                self.reconstruction.object[..., sy, sx] = self.objectPatchUpdate(objectPatch, DELTA)

                # probe update
                self.reconstruction.probe = self.probeUpdate(objectPatch, DELTA)

        end 
    end

end

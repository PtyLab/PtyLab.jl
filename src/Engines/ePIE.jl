export ePIE

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
    createUpdateFunctions(engine::ePIE{T}, objectPatch, probe, DELTA) where T

Return the functions `probeUpdate` and `objectPatchUpdate` which are used
to update the `probe` and `objectPatch`.

We create a lot of memory buffers.
Note, for the `sum` operations, the dimensions over summation is already
specified by the buffer.
"""
function createUpdateFunctions(engine::ePIE{T}, objectPatch, probe, DELTA) where T
    # the let statement is used such that the variables are available in the
    # functions as closures. The actualy operation is not meaningful but
    # instead we want to be sure to have the right shape 
    probeUpdate = let   fracProbe = similar(objectPatch) 
                        newProbe = similar(probe) 
                        fracProbeDELTA = fracProbe .* DELTA
                        sumBufferNewProbe = sum(fracProbeDELTA, dims=(3,4,6)) 
                        abs2objectPatch = real.(objectPatch)
                        sumabs2objectPatch = sum(abs2objectPatch, dims=(3,4,5,6))
        function probeUpdate(objectPatch, probe, DELTA) where T
            abs2objectPatch .= abs2.(objectPatch)
            fracProbe .= conj.(objectPatch) ./ maximum(sum!(sumabs2objectPatch, abs2objectPatch))
            fracProbeDELTA .= fracProbe .* DELTA
            newProbe .= probe .+ engine.betaProbe .* sum!(sumBufferNewProbe, fracProbeDELTA)
            return newProbe
        end
    end
   
    objectPatchUpdate = let   fracObject = similar(probe)
        newObject = similar(objectPatch) 
                        fracObjectDELTA = fracObject .* DELTA
                        sumBufferNewObject = sum(fracObjectDELTA, dims=(3,5,6))
                        abs2probe = real.(probe)
                        sumabs2probe = sum(abs2probe, dims=(3,4,5,6))
        function objectPatchUpdate(objectPatch, probe, DELTA) 
            abs2probe .= abs2.(probe)
            fracObject .= conj.(probe) ./ maximum(sum!(sumabs2probe, abs2probe))
            fracObjectDELTA .= fracObject .* DELTA
            newObject .= objectPatch .+ engine.betaObject .* sum!(sumBufferNewObject, fracObjectDELTA)
            return newObject
        end
    end

    return objectPatchUpdate, probeUpdate
end



"""
    reconstruct(engine::ePIE{T}, params::Params, rec::ReconstructionCPM{T}) where T 

Reconstruct a CPM dataset.
"""
function reconstruct(engine::ePIE{T}, params::Params, rec::ReconstructionCPM{T}) where T 
    # calculate the positions since rec.positions is in fact every time
    # recalculated when called! 
    positions = rec.positions
    # more alias
    object = rec.object
    probe = rec.probe
    Np = rec.Np
    
    intensityProjection!, objectPatchUpdate, probeUpdate, ptychogram, oldProbe, oldObjectPatch, esw, DELTA = 
        _prepareBuffersAndFunctionsPIE(rec, params, engine)
   
    # loop for iterations
    @showprogress for loop in 1:engine.numIterations
        # random order of positions or not
        posList = params.randPositionOrder ? randperm(size(positions, 2)) : (1:size(p, 2))
        
        # critical optimization steps
        _loopUpdatePIE!(posList, positions, Np, ptychogram, 
                        object, oldObjectPatch, oldProbe, probe, esw, DELTA,
                        intensityProjection!, probeUpdate, objectPatchUpdate)
        
        # enforce some constraints like COM
        enforceConstraints!(rec, params)

    end
    return probe, object
end

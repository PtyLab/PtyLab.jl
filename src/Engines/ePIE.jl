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
    reconstruct(engine::ePIE{T}, params::Params, rec::ReconstructionCPM{T}) where T 

Reconstruct a CPM dataset.
"""
function reconstruct(engine::ePIE{T}, params::Params, rec::ReconstructionCPM{T}) where T 
    # calculate the positions since rec.positions is in fact every time
    # recalculated when called! 
    positions = rec.positions::Array{Int64, 2}
    # more alias
    object = rec.object::AbstractArray{Complex{T}, 6}
    probe = rec.probe::AbstractArray{Complex{T}, 6}
    Np = rec.Np
    
    intensityProjection!, objectPatchUpdate, probeUpdate, ptychogram::AbstractArray{T, 3}, oldProbe::typeof(probe), 
        oldObjectPatch::typeof(object), esw::typeof(probe), DELTA::typeof(probe) = 
        _prepareBuffersAndFunctionsPIE(rec, params, engine)
  
    @showprogress for loop in 1:engine.numIterations
        # critical optimization steps
        _loopUpdatePIE!(params.randPositionOrder, positions, Np, ptychogram, 
                        object, oldObjectPatch, oldProbe, probe, esw, DELTA,
                        intensityProjection!, probeUpdate, objectPatchUpdate)
        
        # enforce some constraints like COM
        enforceConstraints!(rec, params)

    end
    return probe, object
end



function f(randPositionOrder, positions, Np, ptychogram, 
                        object, oldObjectPatch, oldProbe, probe, esw, DELTA,
                        intensityProjection!, probeUpdate, objectPatchUpdate, numIterations)
    # loop for iterations
    @showprogress for loop in 1:numIterations
        # critical optimization steps
        _loopUpdatePIE!(randPositionOrder, positions, Np, ptychogram, 
                        object, oldObjectPatch, oldProbe, probe, esw, DELTA,
                        intensityProjection!, probeUpdate, objectPatchUpdate)
        
        # enforce some constraints like COM
        #enforceConstraints!(rec, params)

    end

    return probe, object

end

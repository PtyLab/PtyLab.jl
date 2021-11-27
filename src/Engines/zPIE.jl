export zPIE

"""
    @with_kw mutable struct zPIE{T} <: Engines

`zPIE` is a struct to store important parameters for the `zPIE::Engine`.

"""
@with_kw mutable struct zPIE{T} <: Engines
    betaProbe::T=0.25f0
    betaObject::T=0.25f0
    numIterations::Int=50
    DoF::T = 1
    # gradient step size for axial position correction (typical range [1, 100])
    zPIEgradientStepSize::T = 100  
    zPIEfriction::T = 0.7f0
    focusObject::Bool = True
    zMomentum::T = 0
end



"""
    reconstruct(engine::zPIE{T}, params::Params, rec::ReconstructionCPM{T}) where T 

Reconstruct a CPM dataset.
"""
function reconstruct(engine::zPIE{T}, params::Params, rec::ReconstructionCPM{T}) where T 
    # calculate the positions since rec.positions is in fact every time
    # recalculated when called! 
    positions = rec.positions
    # more alias
    object = rec.object
    probe = rec.probe
    Np = rec.Np
    
    intensityProjection!, objectPatchUpdate, probeUpdate, ptychogram, oldProbe, oldObjectPatch, esw, DELTA = 
        _prepareBuffersAndFunctionsPIE(rec, params, engine)
    
    nlambda_mid = rec.nlambda ÷ 2 + 1

    # loop for iterations
    @showprogress for loop in 1:engine.numIterations

        if loop == 1
            zNew = reconstructio
        else
            dz = range(-rec.DoF, rec.DoF, length=11)
            merit = Vector{T}[]
            
            object_tmp = rec.object[end÷2 - rec.Np÷2 + 1 : end ÷ 2 + n ÷ 2, end÷2 - rec.Np÷2 + 1 : end ÷ 2 + n ÷ 2, nlambda,1,1,1]
            @warn "Fix ASWP prop such that it is outside of loop!"
            aspw_prop = ASPW(object_tmp, rec.zo, rec.spectralDensity[nlambda_mid], rec.Lp)  
            for k = 1:length(dz)
                imProp =  aspw_prop(imgProp)
                aleph = 1f-2
                gradx = imProp .- circshift(imProp, (0, 1))
                grady = imProp .- circshift(imProp, (1, 0))
                merit = push!(merit, sum(sqrt.( abs2.(gradx) + abs2.(grady) .+ aleph)))
            end
            zStep = sum(dz .* merit) / sum(merit);
            eta = 0.7;
            zMomentum = eta * zMomentum + obj.params.zPIEstepSize * zStep;
            zNew = obj.zo + zMomentum;
        end

        rec.z = zNew

        # critical optimization steps
        _loopUpdatePIE!(params.randPositionOrder, positions, Np, ptychogram, 
                        object, oldObjectPatch, oldProbe, probe, esw, DELTA,
                        intensityProjection!, probeUpdate, objectPatchUpdate)
        
        # enforce some constraints like COM
        enforceConstraints!(rec, params)

    end
    return probe, object
end

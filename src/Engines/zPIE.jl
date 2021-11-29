export zPIE

"""
    @with_kw mutable struct zPIE{T} <: Engines

`zPIE` is a struct to store important parameters for the `zPIE::Engine`.

"""
@with_kw mutable struct zPIE{T} <: Engines
    betaProbe::T=0.25f0
    betaObject::T=0.25f0
    numIterations::Int=50
    DoF::T = 1f0
    # gradient step size for axial position correction (typical range [1, 100])
    zPIEgradientStepSize::T = 100f0 
    zPIEfriction::T = 0.7f0
    focusObject::Bool = true
    zMomentum::T = 0f0
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

    zMomentum = zero(T)
    # loop for iterations
    aspw_prop = ASPW(esw[:, :, 1,1,1,1], rec.zo, rec.spectralDensity[nlambda_mid], rec.Lp)  
    
    @showprogress for loop in 1:engine.numIterations

        if loop == 1
            zNew = rec.zo 
        else
            dz = range(-rec.DoF, rec.DoF, length=11)
            # vector of Float32, Float64, ... 
            merit = T[]
            
            object_tmp = rec.object[end÷2 - rec.Np÷2 + 1 : end ÷ 2 + rec.Np ÷ 2, 
                                    end÷2 - rec.Np÷2 + 1 : end ÷ 2 + rec.Np ÷ 2, rec.nlambda,1,1,1]
            for k = 1:length(dz)
                imProp =  aspw_prop(copy(object_tmp), ASPW_kernel(object_tmp, dz[k], rec.spectralDensity[nlambda_mid], rec.Lp))
                aleph = 1f-2
                gradx = imProp .- circshift(imProp, (0, 1))
                grady = imProp .- circshift(imProp, (1, 0))
                value_merit = sum(sqrt.( abs2.(gradx) + abs2.(grady) .+ aleph))
                #@show value_merit
                merit = push!(merit, value_merit)
            end
            zStep = sum(dz .* merit) / sum(merit);
            eta = 0.7;
            zMomentum = eta * zMomentum + engine.zPIEgradientStepSize * zStep;
            zNew = rec.zo + zMomentum;
            @show zNew
        end

        rec.zo = zNew

        # critical optimization steps
        _loopUpdatePIE!(params.randPositionOrder, positions, Np, ptychogram, 
                        object, oldObjectPatch, oldProbe, probe, esw, DELTA,
                        intensityProjection!, probeUpdate, objectPatchUpdate)
        
        # enforce some constraints like COM
        enforceConstraints!(rec, params)

    end
    return probe, object
end

"""
    reconstruct(engine::Union{mPIE{T}, ePIE{T}}, params::Params, rec::ReconstructionCPM{T}) where T 

Reconstruct a CPM dataset.


"""
function reconstruct(engine::Union{mPIE{T}, ePIE{T}}, params::Params, rec::ReconstructionCPM{T}) where T 
    # calculate the positions since rec.positions is in fact every time
    # recalculated when called! 
    positions = rec.positions::Array{Int64, 2}
    # more alias
    object = rec.object::AbstractArray{Complex{T}, 6}
    probe = rec.probe::AbstractArray{Complex{T}, 6}
    Np = rec.Np
    
    intensityProjection!, updateFunctions, ptychogram::AbstractArray{T, 3}, oldProbe::typeof(probe), 
        oldObjectPatch::typeof(object), esw::typeof(probe), DELTA::typeof(probe) = 
        _prepareBuffersAndFunctionsPIE(rec, params, engine)
  
    @showprogress for loop in 1:engine.numIterations
        # critical optimization steps
        _loopUpdatePIE!(engine, params.randPositionOrder, positions, Np, ptychogram, 
                        object, oldObjectPatch, oldProbe, probe, esw, DELTA,
                        intensityProjection!, updateFunctions)
        
        # enforce some constraints like COM
        # causes some allocations since circshift causes them 
        # for higher dimensional arrays
        enforceConstraints!(rec, params)

    end
    return probe, object
end





"""
    _loopUpdatePIE!(<args>)

Internal function which serves as the inner for loop of the PIE based algorithms
like `ePIE` or `zPIE`.
    
"""
function _loopUpdatePIE!(engine, randPositionOrder, positions, Np, ptychogram, 
                         object, oldObjectPatch, oldProbe, probe, esw, DELTA,
                         intensityProjection!, updateFunctions)


    # random order of positions or not
    posList = randPositionOrder ? randperm(size(positions, 2)) : (1:size(positions, 2))

    @assert size(positions, 1) == 2 "encoder has shape $(size(encoder)) but needs to have tranposed shape!"
    for positionIndex in posList 
        row, col = view(positions, :, positionIndex) 

        sy = row:(row + Np - 1)
        sx = col:(col + Np - 1)
        # no copy necessary -> because esw immediately calculated
        objectPatch = view(object, sy, sx, ..)
        
        # save old state since we need that later for probeUpdate and objectPatchUpdate
        oldObjectPatch .= objectPatch
        oldProbe .= probe
        
        # exit surface wave,
        esw .= objectPatch .* probe
        
        # already store esw in DELTA, since intensityProjection is going to change esw
        DELTA .= -1 .* esw
        eswUpdate = intensityProjection!(esw, view(ptychogram, :, :, positionIndex))
        
         # difference term
        DELTA .+= eswUpdate
        
        # update newProbe und newObjectPatch
        newObjectPatch = updateFunctions.objectPatchUpdate(oldObjectPatch, oldProbe, DELTA)
        newProbe = updateFunctions.probeUpdate(oldObjectPatch, oldProbe, DELTA)
        probe .= newProbe
        object[sy, sx, ..] .= newObjectPatch

        
        if typeof(engine) == mPIE
            if rand(1) > 0.95
                updateFunctions.objectMomentumUpdate(object)
                updateFunctions.probeMomentumUpdate(probe)
            end
        end
    end 

    return nothing
end

"""
    _prepareBuffersAndFunctionsPIE(rec, params, engine)

Prepate buffers and functions needed for PIE based reconstructions
"""
function _prepareBuffersAndFunctionsPIE(rec, params, engine::Union{ePIE, mPIE})
    # create intensityProjection function 
    intensityProjection! = IntensityProjection(rec, params)

    # get two function which maybe shift depending on the flag             \
    maybe_fftshift, maybe_ifftshift = getMaybeFftshifts(! params.fftshiftFlag) 

    # ptychogram = rec.ptychogram
    ptychogram = maybe_ifftshift(rec.ptychogram)

    # if that flag is true, swap first and second dimension (transpose)
    if params.transposePtychogram
        ptychogram = permutedims(ptychogram, (2,1,3))
    end

    # those are buffers
    oldProbe = copy(rec.probe)
    oldObjectPatch = rec.object[1:rec.Np, 1:rec.Np, ..]
    # just a buffer with correct shape
    esw = rec.object[1:rec.Np, 1:rec.Np, ..] .* rec.probe
    DELTA = similar(esw)
    # performant update functions for probe and patch without memory usage.
    updateFunctions = createUpdateFunctions(engine, oldObjectPatch, oldProbe, rec.object, DELTA)

    return (intensityProjection!, updateFunctions,
            ptychogram, oldProbe, oldObjectPatch, esw, DELTA)
end




export reconstruct

abstract type Engines end

"""
    centerOfMassStabilizationOffset(probe)

Calculate the center of mass in the first 2 dimensions 
of a N dimensional array `probe`.
Return a tuple indicating the center.
"""
function centerOfMassStabilizationOffset(probe)
    # calculate center of mass 
    P2 = sum(abs2, probe, dims=(3,4,5,6));
    # total weight of the probe
    denom = sum(P2)
    # weight each pixel with a distance
    xc = round(Int, sum((1:size(probe, 2)) .* P2) ./ denom);
    yc = round(Int, sum((1:size(probe, 1))' .* P2) ./ denom);
    return (xc, yc) 
end



function enforceConstraints!(rec::ReconstructionCPM{T}, params)  where T
    # enforce center of mass of probe
    if params.comStabilizationSwitch
        # fourier center (but not important where exactly)
        center = (size(rec.probe, 1), size(rec.probe, 2)) .÷ 2 .+ 1
        offset = center .- centerOfMassStabilizationOffset(rec.probe)
        # shift in direction to the center
        rec.probe .= circshift(rec.probe, offset);
    end
end



"""
    IntensityProjection(rec::ReconstructionCPM{T}, params::Params)

Returns a function which does the `intensityProjection!`.
The returned function will mutate the input `esw`.

 # Examples
```julia
julia> intensityProjection! = IntensityProjection(rec, params)

julia> eswUpdate = intensityProjection!(esw, Imeasured)
```
"""
function IntensityProjection(rec::ReconstructionCPM{T}, params::Params) where T
    # create an efficient propagator function
    esw_temp = rec.object[1:rec.Np, 1:rec.Np, ..] .* rec.probe
    object2detector, detector2object = params.propagatorType(esw_temp)
      
    
    @warn "gimmel is currently estimated as `100 * eps($T)`"
    gimmel = 100 * eps(T)
    f! = let    intensityConstraint = params.intensityConstraint
                object2detector = object2detector
                detector2object = detector2object
                abs2_buffer = real(similar(esw_temp))
                sum_buffer = real(similar(esw_temp, (size(esw_temp, 1), size(esw_temp, 2), 1, 1, 1, 1)))
                frac_buffer = similar(sum_buffer) 
		# TODO, check if better solution exists
    		Iestimated = similar(real.(rec.probe), (size(rec.probe, 1), size(rec.probe, 2)))
        function f(esw, Imeasured)
            ESW = object2detector(esw)
            Iestimated = let
                if intensityConstraint === IntensityConstraintStandard
                    # sum over the last three channels.
                    map!(abs2, abs2_buffer, ESW)
                    view(sum!(sum_buffer, abs2_buffer), :, :, 1,1,1,1)
                else
                    error("Unknown intensityConstraint $intensityConstraint")
                end
            end

            frac = let 
                if intensityConstraint === IntensityConstraintStandard 
                     frac_buffer .= sqrt.(Imeasured ./ (Iestimated .+ gimmel))
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
    _loopUpdatePIE!(<args>)

Internal function which serves as the inner for loop of the PIE based algorithms
like `ePIE` or `zPIE`.
    
"""
function _loopUpdatePIE!(posList, positions, Np, ptychogram, 
                         object, oldObjectPatch, oldProbe, probe, esw, DELTA,
                         intensityProjection!, probeUpdate, objectPatchUpdate)
    for positionIndex in posList 
        row, col = positions[:, positionIndex] 
            
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
        newProbe = probeUpdate(oldObjectPatch, oldProbe, DELTA)
        newObjectPatch = objectPatchUpdate(oldObjectPatch, oldProbe, DELTA)
        probe .= newProbe
        object[sy, sx, ..] .= newObjectPatch
    end 
end

"""
    _prepareBuffersAndFunctionsPIE(rec, params, engine)

Prepate buffers and functions needed for PIE based reconstructions
"""
function _prepareBuffersAndFunctionsPIE(rec, params, engine)
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
    objectPatchUpdate, probeUpdate = createUpdateFunctions(engine, oldObjectPatch, oldProbe, DELTA)

    return (intensityProjection!, objectPatchUpdate, probeUpdate, 
            ptychogram, oldProbe, oldObjectPatch, esw, DELTA)
end

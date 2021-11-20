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
    object2detector, detector2object = params.propagatorType(esw_temp)
       
    
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

export reconstruct


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
    esw_temp = view(rec.object, 1:rec.Np, 1:rec.Np, ..) .* rec.probe
    object2detector, detector2object = params.propagatorType(esw_temp)
      

    # if those two variables are set, we are only allowed to apply the updates
    # on known sensor values
    multiSensor = (!isnothing(rec.Nd1)) && (!isnothing(rec.Nd2))

    mask = let
        if multiSensor
            arr = similar(rec.ptychogram, Bool, rec.Nd1, rec.Nd2, length(rec.posDetectors))
            fill!(arr, 1)
            mask = assembleMultiSensorPtychogram(arr, rec.Nd, rec.Nd1, 
                                                 rec.Nd2, rec.posDetectors, rec.Ld, rec.dxd)[:, :, 1]
            mask = ifftshift(mask)
        else
            nothing
        end
    end


    @warn "gimmel is currently estimated as `eps($T)`"
    f! = let    intensityConstraint = params.intensityConstraint
                gimmel = eps(T)
                abs2_buffer = similar(esw_temp, real(eltype(esw_temp)))
                sum_buffer = similar(esw_temp, real(eltype(esw_temp)), (size(esw_temp, 1), size(esw_temp, 2), 1, 1, 1, 1))
                frac_buffer = similar(sum_buffer) 
        function f(esw, Imeasured)
            ESW = object2detector(esw)
            Iestimated = let
                if typeof(intensityConstraint) === IntensityConstraintStandard
                    # sum over the last three channels.
                    map!(abs2, abs2_buffer, ESW)
                    view(sum!(sum_buffer, abs2_buffer), :, :, 1,1,1,1)
                else
                    error("Unknown intensityConstraint $intensityConstraint")
                end
            end

            frac = let 
                if typeof(intensityConstraint) === IntensityConstraintStandard 
                     frac_buffer .= sqrt.(Imeasured ./ (Iestimated .+ gimmel))
                else
                    error("Unknown intensityConstraint")
                end
            end

            # update ESW and if there is as mask (for multiple sensors)
            if isnothing(mask)
                ESW .*= frac 
            else
                # frac is only multiplied to those parts of the array
                # where we have data available to compare with
                # the other parts are not explicitly updated
                # but still they are changed implicitly with forward and
                # backward propagations

                # this line causes allocations
                # https://github.com/JuliaLang/julia/issues/43442
                ESW[mask, ..] .*= view(frac, mask, ..)
            end

            # back to detector, memory free due to plan_fft!
            eswUpdate = detector2object(ESW)
            return eswUpdate 
        end
    end

    return f! 
end

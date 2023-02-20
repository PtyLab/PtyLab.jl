export simshow
export show_grid

"""
    simshow(arr; set_zero=false, set_one=false, γ=1, cmap=:gray)
Displays a real valued array.
Works within Jupyter and Pluto.
# Keyword args
The transforms are applied in that order.
* `set_zero=false` subtracts the minimum to set minimum to 0
* `set_one=true` divides by the maximum to set maximum to 1
* `γ` applies a gamma correction
* `cmap=:gray` applies a colormap provided by ColorSchemes.jl. If `cmap=:gray` simply `Colors.Gray` is used
    and with different colormaps the result is an `Colors.RGB` element type.  
    You can try `:jet`, `:deep`, `thermal` or different ones by reading the catalogue of ColorSchemes.jl
"""
function simshow(arr::AbstractArray{T};
                 set_zero=false, set_one=true,
                 γ = one(T),
                 cmap=:gray) where {T<:Real}

    arr = set_zero ? arr .- minimum(arr) : arr

    if set_one
        m = maximum(arr)
        if !iszero(m)
            arr = arr ./ maximum(arr)
        end
    end

    if !isone(γ)
        arr = arr .^ γ
    end

    if cmap == :gray
        Gray.(arr)
    else
        get(colorschemes[cmap], arr)
    end
end


"""
    simshow(arr; γ=1)
Displays a complex array. Color encodes phase, brightness encodes magnitude.
Works within Jupyter and Pluto.
# Keyword args
* `γ` applies a gamma correction to the magnitude
"""
function simshow(arr::AbstractArray{T};
                 γ=one(T), kwargs...) where (T<:Complex)

    Tr = real(T)
    # scale abs to 1
    absarr = abs.(arr)
    absarr ./= maximum(absarr)

    if !isone(γ)
        absarr .= absarr .^ γ
    end

    angarr = angle.(arr) ./ Tr(2pi) * Tr(360)

    HSV.(angarr, one(Tr), absarr)
end


"""
    show_grid(grr::PtyLab.GridRegularRand; only_points=false, thickness=2)
Displays the `PtyLab.GridRegularRand`.
`only_points=true` displays anchor points otherwise lines.
`thickness` is the thickness of the drawing
"""
function show_grid(grr::PtyLab.GridRegularRand; only_points=false, thickness=1)
    img = zeros(grr.grid_size)
    # draw lines for the grid
    for t in grr.tiles
        if only_points
            img[t.i₁:t.i₁+thickness, t.j₁:t.j₁+thickness] .= 1
        else
            img[t.i₁:t.i₂, t.j₁:t.j₁+thickness-1] .= 1
            img[t.i₁:t.i₂, t.j₂:t.j₂+thickness-1] .= 1
            img[t.i₁:t.i₁+thickness-1, t.j₁:t.j₂] .= 1
            img[t.i₂:t.i₂+thickness-1, t.j₁:t.j₂] .= 1
        end
    end
    simshow(img)
end

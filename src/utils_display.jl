export complex_show, gray_show, show_grid 

"""
    complex_show(arr)

Displays a complex array. Color encodes phase, brightness encodes magnitude.
Works within Jupyter and Pluto.
"""
function complex_show(cpx::AbstractArray{T, N}) where {T<:Complex, N}
    Tr = real(T)
	ac = abs.(cpx)
    HSV.(angle.(cpx)./Tr(2pi)*256, ones(Tr,size(cpx)), ac./maximum(ac))
end


"""
    gray_show(arr; set_one=false, set_zero=false)

Displays a real gray color array. Brightness encodes magnitude.
Works within Jupyter and Pluto.

## Keyword args
* `set_one=false` divides by the maximum to set maximum to 1
* `set_zero=false` subtracts the minimum to set minimum to 1
"""
function gray_show(arr; set_one=true, set_zero=false)
    arr = set_zero ? arr .- minimum(arr) : arr
    arr = set_one ? arr ./ maximum(arr) : arr
    Gray.(arr)
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
    gray_show(img) 
end

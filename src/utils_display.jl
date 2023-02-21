export show_grid

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

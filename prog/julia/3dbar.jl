using GLMakie
using GeometryBasics

function Bar3D(x, y, z)
    cbarPal = :Spectral_11
    fig = Figure(resolution = (800, 600))
    ax  = Axis3(fig[1, 1]; aspect = (1,1,1), elevation = π/6, perspectiveness = 0.5)

    n = length(x)
    δx = (x[2] - x[1]) / 2
    δy = (y[2] - y[1]) / 2
    rectMesh = Rect3f(Vec3f(-0.5, -0.5, 0), Vec3f(1, 1, 1))
    
    meshscatter!(ax, x, y, 0 * z;
        marker = rectMesh,
        color = z[:],
        markersize = Vec3f.(2δx, 2δy, z[:]),
        # colormap = cbarPal,
    )
    limits!(ax, -1, n + 1, -1, n + 1, -2, 2)
    fig
end

x = y = 1:5
Bar3D(x, y, randn(5,5))
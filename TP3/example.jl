using Makie
using PlotlyJS

xs = LinRange(-1, 1, 100)
ys = LinRange(-1, 1, 100)
zs = [x^2+y^2 for x in xs, y in ys]

s = contour(xs, ys, zs, transformation=(:yz, -1), levels=10, limits=FRect((-2,-2,-2), (4,4,4)))
contour!(s, xs, ys, zs, transformation=(:yz,  0), levels=10)
contour!(s, xs, ys, zs, transformation=(:yz,  1), levels=10)
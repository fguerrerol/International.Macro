using Makie
using JLD2
using AbstractPlotting

A = load_object("SOE_cc.jdl2")
# Parent scence
scene, layout = layoutscene()

# Scene for main plot
lscene = layout[1, 1] = LScene(scene, camera=cam3d!, raw=false)
s = lscene.scene

# Main plot

style = (levels=10, linewidth=2)
S = contour(A.zgrid,A.ξgrid,A.Y[1,:,:], transformation=(:yz, -1); style..., limits=FRect((-2,-2,-2), (4,4,4)))
contour!(S,A.zgrid,A.ξgrid,A.Y[2,:,:], transformation=(:yz,  0); style...)
contour!(S, A.zgrid,A.ξgrid,A.Y[3,:,:], transformation=(:yz,  1); style...)

# Axis names
xlabel!(s, "X [mm]")
ylabel!(s, "Y [mm]")
zlabel!(s, "Z⋅10^3 [mm]")
s[Axis][:names][:rotation] = (qrotation(Vec3(0,0,1), π),
                              qrotation(Vec3(0,0,1), π/2),
                              qrotation(Vec3(0,0,1), -π/2) * qrotation(Vec3(0,1,0), -π/2))
s[Axis][:names][:align] = ((:center, :right), (:center, :right), (:center, :left))

# Axis ticks
xtickrotation!(s, 1π)   # Giving an Irrational fails (probably a bug)
ytickrotation!(s, π/2)

# Colorbar
layout[1, 2] = LColorbar(scene, s[end], width=30)

scene

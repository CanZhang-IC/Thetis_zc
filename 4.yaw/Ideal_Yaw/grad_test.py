from thetis import *

lx = 40e3
ly = 2e3
nx = 25
ny = 2
mesh2d = RectangleMesh(nx, ny, lx, ly)

x = SpatialCoordinate(mesh2d)

n = as_vector((sin(x[0]/180*pi),cos(x[0]/180*pi)))
for i in [0,1]:
    for j in [0,1]:
        print(i,j)
        a = assemble(grad(n)[i,j]*dx)
        print(a)
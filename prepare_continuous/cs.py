from thetis import *
from firedrake_adjoint import *
import myboundary



mesh2d = Mesh('./mesh/mesh.msh')
#timestepping options
dt = 30*60 # reduce this if solver does not converge
t = 215*30*60


P1 = FunctionSpace(mesh2d, "CG", 1)
VP = VectorFunctionSpace(mesh2d, "CG", 1)

tidal_v = Function(VP)
v = myboundary.set_velocity_field(Function(VP), t,dt)
tidal_v.project(v)
File('v.pvd').write(tidal_v) 


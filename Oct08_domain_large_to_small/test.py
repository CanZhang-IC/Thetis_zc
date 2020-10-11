from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
op2.init(log_level=INFO)
from os.path import join
import detectors
import tidal_forcing
import utm
import numpy
import scipy.interpolate
import myboundary

ouput_dir = 'outputs/see-ele'


mesh2d = Mesh('./mesh/mesh.msh')
P1 = FunctionSpace(mesh2d, "CG", 1)
tidal_elev = Function(P1)
elev = myboundary.set_tidal_field(Function(P1), 2, 1)
tidal_elev.project(elev)
File('ele.pvd').write(tidal_elev)

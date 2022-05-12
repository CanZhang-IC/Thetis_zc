"""
An ideal case for yaw angle optimisation
"""
# to enable a gradient-based optimisation using the adjoint to compute
# gradients, we need to import from thetis_adjoint instead of thetis. This
# ensure all firedrake operations in the Thetis model are annotated
# automatically, in such a way that we can rerun the model with different input
# parameters, and also derive the adjoint-based gradient of a specified input
# (the functional) with respect to a specified input (the control)
from ast import Constant
from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
op2.init(log_level=INFO)
import sys
sys.path.append('..')
import os
import time

import yagmail
import matplotlib.pyplot as plt
import numpy as np
import h5py



### set up the Thetis solver obj as usual ##
mesh2d = Mesh('./mesh.msh')


#set viscosity bumps at in-flow boundaries.
P1_2d = FunctionSpace(mesh2d, 'CG', 1)
V = FunctionSpace(mesh2d, "CG", 1)
DG_2d = FunctionSpace(mesh2d, "DG", 1)
vector_dg = VectorFunctionSpace(mesh2d, "DG", 1)

for i in range(41):
    while len(str(i)) < 5:
        i = '0'+str(i)
    # initialise velocity and elevation
    chk = DumbCheckpoint("../../sediment_hydro/outputs-hydro/pure-hydro-ss50/hdf5/Elevation2d_"+i, mode=FILE_READ)
    elev1 = Function(DG_2d, name="elev_2d")
    chk.load(elev1)
    chk.close()

    chk = DumbCheckpoint('../../sediment_hydro/outputs-hydro/pure-hydro-ss50/hdf5/Velocity2d_'+i, mode=FILE_READ)
    uv1 = Function(vector_dg, name="uv_2d")
    chk.load(uv1)
    chk.close()

    # initialise velocity and elevation
    chk = DumbCheckpoint("../../sediment_hydro/outputs-hydro/pure-hydro-nosediment/hdf5/Elevation2d_"+i, mode=FILE_READ)
    elev2 = Function(DG_2d, name="elev_2d")
    chk.load(elev2)
    chk.close()

    chk = DumbCheckpoint('../../sediment_hydro/outputs-hydro/pure-hydro-nosediment/hdf5/Velocity2d_'+i, mode=FILE_READ)
    uv2 = Function(vector_dg, name="uv_2d")
    chk.load(uv2)
    chk.close()

    elev3 = Function(DG_2d,name='elevation').project(elev1-elev2)
    File('../../sediment_hydro/diff/elev_diff.pvd',mode='a').write(elev3)

    uv3 = Function(vector_dg,name='uvation').project(uv1-uv2)
    File('../../sediment_hydro/diff/uv_diff.pvd',mode='a').write(uv3)
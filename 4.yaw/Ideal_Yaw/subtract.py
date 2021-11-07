from thetis import*
from firedrake import*
import numpy as np
import matplotlib.pyplot as plt


mesh2d = Mesh('headland2.msh')

P1_2d = FunctionSpace(mesh2d, 'DG', 1)
V = VectorFunctionSpace(mesh2d,'DG',1)

checkpoint_file1 = DumbCheckpoint('./outputs_60/hdf5/Velocity2d_00039', mode = FILE_READ)
uv_init_1 = Function(V,name='uv_2d')
checkpoint_file1.load(uv_init_1)

checkpoint_file2 = DumbCheckpoint('./outputs_0/hdf5/Velocity2d_00039', mode = FILE_READ)
uv_init_2 = Function(V,name='uv_2d')
checkpoint_file2.load(uv_init_2)

subtract_u = Function(V,name='uv_2d').interpolate(uv_init_2 - uv_init_1)

File('subtract_u3.pvd').write(subtract_u)
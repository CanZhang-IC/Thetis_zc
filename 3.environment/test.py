from thetis import*
from firedrake import*
import numpy as np
import matplotlib.pyplot as plt


tide_numbers = [504,516]

for tide_number in range(492,519):


	mesh2d = Mesh('../mesh/mesh.msh')

	P1_2d = FunctionSpace(mesh2d, 'DG', 1)
	V = VectorFunctionSpace(mesh2d,'DG',1)

	checkpoint_file1 = DumbCheckpoint('../../outputs/2.economy/discrete/flood_ebb-forward/cable-BE-1/hdf5/Velocity2d_00'+str(tide_number), mode = FILE_READ)
	uv_init_1 = Function(V,name='uv_2d')
	checkpoint_file1.load(uv_init_1)

	checkpoint_file2 = DumbCheckpoint('../../outputs/0.validation/discrete-4cores/hdf5/Velocity2d_00'+str(tide_number), mode = FILE_READ)
	uv_init_2 = Function(V,name='uv_2d')
	checkpoint_file2.load(uv_init_2)

	uv3 = Function(V,name='uv_2d').interpolate(-uv_init_1+uv_init_2)
	File('diff/velocity_diff.pvd', mode = 'a').write(uv3)






		




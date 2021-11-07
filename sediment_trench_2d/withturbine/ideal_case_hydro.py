"""
An ideal case for yaw angle optimisation
"""
# to enable a gradient-based optimisation using the adjoint to compute
# gradients, we need to import from thetis_adjoint instead of thetis. This
# ensure all firedrake operations in the Thetis model are annotated
# automatically, in such a way that we can rerun the model with different input
# parameters, and also derive the adjoint-based gradient of a specified input
# (the functional) with respect to a specified input (the control)
from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
op2.init(log_level=INFO)
import sys
sys.path.append('..')
import os
import numpy
import time
import h5py
import yagmail

t_start = time.time()

H = 40
distance = 10
speed = 2
output_dir = '../../../outputs/sediment/test/hydro'

### set up the Thetis solver obj as usual ##
mesh2d = Mesh('../../prepare_ideal_meshes/conference_mesh2_with_effected_area.msh')

dt = 30
t_end = dt*100
# export interval in seconds
t_export = np.round(t_end/40, 0)

#set viscosity bumps at in-flow boundaries.
P1_2d = FunctionSpace(mesh2d, 'CG', 1)
x = SpatialCoordinate(mesh2d)
h_viscosity = Constant(1e-6)

# create solver and set options
solver_obj = solver2d.FlowSolver2d(mesh2d, Constant(H))
options = solver_obj.options
options.timestep = dt
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.output_directory = output_dir
options.check_volume_conservation_2d = True
options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']
options.element_family = 'dg-dg'
options.timestepper_type = 'CrankNicolson'
options.timestepper_options.implicitness_theta = 1.0
options.timestepper_options.use_semi_implicit_linearization = True
# using direct solver as PressurePicard does not work with dolfin-adjoint (due to .split() not being annotated correctly)
options.timestepper_options.solver_parameters = {'snes_monitor': None,
                                                'snes_rtol': 1e-5,
                                                'ksp_type': 'preonly',
                                                'pc_type': 'lu',
                                                'pc_factor_mat_solver_type': 'mumps',
                                                'mat_type': 'aij'
                                                }
options.horizontal_viscosity = h_viscosity
# options.quadratic_drag_coefficient = Constant(0.0025)

# define parameters
average_size = 160*(10**(-6))
ksp = Constant(3*average_size)
# using nikuradse friction
options.nikuradse_bed_roughness = ksp

options.norm_smoother = Constant(0.1)

# assign boundary conditions
left_tag = 1
right_tag = 2
coasts_tag = 3
tidal_elev = Function(get_functionspace(mesh2d, "CG", 1), name='tidal_elev').assign(0.0)
tidal_vel = Function(P1_2d).assign(0.0)

solver_obj.bnd_functions['shallow_water'] = {
    right_tag: {'uv':as_vector((-Constant(speed), 0.0))},
    left_tag: {'elev': tidal_elev},
}

# def update_forcings(t):
#     print_output("Updating tidal elevation at t = {}".format(t))
#     tidal_vel.project(-Constant(speed))

#set initial condition
solver_obj.assign_initial_conditions(uv=as_vector((-Constant(speed), 0.0)), elev=tidal_elev)

# start computer forward model
solver_obj.iterate()#update_forcings=update_forcings)

uv, elev = solver_obj.fields.solution_2d.split()

checkpoint_dir = "../../../outputs/sediment/hydrodynamics_trench"

# if not os.path.exists(checkpoint_dir):
#     os.makedirs(checkpoint_dir)
chk = DumbCheckpoint(checkpoint_dir + "/velocity", mode=FILE_CREATE)
chk.store(uv, name="velocity")
chk.close()
chk = DumbCheckpoint(checkpoint_dir + "/elevation", mode=FILE_CREATE)
chk.store(elev, name="elevation")
chk.close()

t_end = time.time()
print('time cost: {0:.2f}min'.format((t_end - t_start)/60))




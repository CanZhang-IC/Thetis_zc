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

t_start = time.time()

H = 40
output_dir = '../../../outputs/4.yaw/Conference2/restart'

### set up the Thetis solver obj as usual ###
mesh2d = Mesh('../../prepare_ideal_meshes/headland2.msh')

tidal_amplitude = 5.
tidal_period = 12.42*60*60
timestep = 300
t_export = 600
t_end =  tidal_period * 3 

#set viscosity bumps at in-flow boundaries.
P1_2d = FunctionSpace(mesh2d, 'CG', 1)
x = SpatialCoordinate(mesh2d)
v_b = 50
v_inner = 0.1
v_length = 200
lx = 2000
h_viscosity = Function(P1_2d).interpolate(conditional(le(x[0], v_length), v_b+v_inner-x[0]*v_b/v_length, conditional(ge(x[0],lx-v_length),(x[0]-(lx-v_length))*v_b/v_length+v_inner,v_inner)))
File(output_dir+'/viscosity.pvd').write(h_viscosity)

# create solver and set options
solver_obj = solver2d.FlowSolver2d(mesh2d, Constant(H))
options = solver_obj.options
options.cfl_2d = 1.0
options.use_nonlinear_equations = True
options.timestep = timestep
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.output_directory = output_dir
options.check_volume_conservation_2d = True
options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']
options.element_family = 'dg-dg'
options.timestepper_type = 'CrankNicolson'
options.timestepper_options.implicitness_theta = 1
# using direct solver as PressurePicard does not work with dolfin-adjoint (due to .split() not being annotated correctly)
options.timestepper_options.solver_parameters = {'snes_monitor': None,
                                                 'snes_rtol': 1e-3,
                                                 'ksp_type': 'preonly',
                                                 'pc_type': 'lu',
                                                 'pc_factor_mat_solver_type': 'mumps',
                                                 'mat_type': 'aij'
                                                 }
options.horizontal_viscosity =h_viscosity
options.quadratic_drag_coefficient = Constant(0.0025)

# Boundary conditions - Steady state case
tidal_elev = Function(P1_2d).assign(0.0)
solver_obj.bnd_functions['shallow_water'] = {1: {'elev': tidal_elev},
                                             2: {'elev':tidal_elev}}

# a function to update the tidal_elev bc value every timestep
x = SpatialCoordinate(mesh2d)
g = 9.81
omega = 2 * pi / tidal_period

def update_forcings(t):
    print_output("Updating tidal elevation at t = {}".format(t))
    t = Constant(t)
    tidal_elev.project(tidal_amplitude*sin(omega*t + omega/pow(g*H, 0.5)*x[0]))

#set initial condition
solver_obj.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=tidal_elev)

locations, names = [(1000,300)], ['centre points']
cb = DetectorsCallback(solver_obj, locations, ['elev_2d', 'uv_2d'], name='detectors',detector_names=names)
solver_obj.add_callback(cb,'timestep')

# start computer forward model
solver_obj.iterate(update_forcings=update_forcings)

t_end = time.time()
print('time cost: {0:.2f}min'.format((t_end - t_start)/60))


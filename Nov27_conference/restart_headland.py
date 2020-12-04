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
output_dir = '../../outputs/conference/restart_300_dgdg'

### set up the Thetis solver obj as usual ###
mesh2d = Mesh('../prepare_ideal_meshes/headland2.msh')

tidal_amplitude = 5.
tidal_period = 12.42*60*60
timestep = 300
t_export = 2 * timestep
t_end =  tidal_period + 60*t_export 


#set viscosity bumps at in-flow boundaries.
P1_2d = FunctionSpace(mesh2d, 'CG', 1)
x = SpatialCoordinate(mesh2d)
h_viscosity = Function(P1_2d).interpolate(conditional(le(x[0], 50), 51-x[0], conditional(ge(x[0],1950),x[0]-1949,1)))

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
options.element_family = 'dg-cg'
options.timestepper_type = 'CrankNicolson'
options.timestepper_options.implicitness_theta = 1
# using direct solver as PressurePicard does not work with dolfin-adjoint (due to .split() not being annotated correctly)
options.timestepper_options.solver_parameters = {'snes_monitor': None,
                                                 'snes_rtol': 1e-9,
                                                 'ksp_type': 'preonly',
                                                 'pc_type': 'lu',
                                                 'pc_factor_mat_solver_type': 'mumps',
                                                 'mat_type': 'aij'
                                                 }
options.horizontal_viscosity =h_viscosity
options.quadratic_drag_coefficient = Constant(0.0025)
# options.manning_drag_coefficient = Constant(0.02)

# assign boundary conditions
left_tag = 1
right_tag = 2
coasts_tag = 3
tidal_elev = Function(get_functionspace(mesh2d, "CG", 1), name='tidal_elev')
tidal_elev_bc = {'elev': tidal_elev}

# noslip currently doesn't work (vector Constants are broken in firedrake_adjoint)
freeslip_bc = {'un': Constant(0.0)}
solver_obj.bnd_functions['shallow_water'] = {
    left_tag: tidal_elev_bc,
    right_tag: tidal_elev_bc,
    coasts_tag: freeslip_bc
}

# a function to update the tidal_elev bc value every timestep
x = SpatialCoordinate(mesh2d)
g = 9.81
omega = 2 * pi / tidal_period

def update_forcings(t):
    print_output("Updating tidal elevation at t = {}".format(t))
    tidal_elev.project(tidal_amplitude*sin(omega*t + omega/Constant(pow(g*H, 0.5))*x[0]))

#set initial condition
solver_obj.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=tidal_elev)
# solver_obj.load_state(60, '../../outputs/conference/restart_300000')

locations, names = [(1000,300)], ['centre points']
cb = DetectorsCallback(solver_obj, locations, ['elev_2d', 'uv_2d'], name='detectors',detector_names=names)
solver_obj.add_callback(cb,'timestep')

# start computer forward model
solver_obj.iterate(update_forcings=update_forcings)

t_end = time.time()
print('time cost: {0:.2f}min'.format((t_end - t_start)/60))


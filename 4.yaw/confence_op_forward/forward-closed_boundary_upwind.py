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
output_dir = '../../../outputs/4.yaw/Yaw_Ideal/op-conference_mesh2-5_40/f30-cos00/forward-staggered-both-op'
#Comment for testing forward model

### set up the Thetis solver obj as usual ##
mesh2d = Mesh('../../prepare_ideal_meshes/conference_mesh2.msh')

tidal_amplitude = 5.
tidal_period = 12.42*60*60
timestep = 60
t_export = 2 * timestep
t_end = 20*t_export #12000

#set viscosity bumps at in-flow boundaries.
P1_2d = FunctionSpace(mesh2d, 'CG', 1)
x = SpatialCoordinate(mesh2d)
h_viscosity = Constant(1e-3)

# create solver and set options
solver_obj = solver2d.FlowSolver2d(mesh2d, Constant(H))
options = solver_obj.options
options.timestep = timestep
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
options.quadratic_drag_coefficient = Constant(0.0025)

# assign boundary conditions
left_tag = 1
right_tag = 2
coasts_tag = 3
tidal_elev = Function(get_functionspace(mesh2d, "CG", 1), name='tidal_elev')

tidal_vel = Function(P1_2d).assign(0.0)

solver_obj.bnd_functions['shallow_water'] = {
    right_tag: {'un':tidal_vel},
    left_tag: {'elev': tidal_elev},
}

# a function to update the tidal_elev bc value every timestep
x = SpatialCoordinate(mesh2d)

# Initialise Discrete turbine farm characteristics
farm_options = DiscreteTidalTurbineFarmOptions()
farm_options.turbine_type = 'constant'
farm_options.turbine_options.thrust_coefficient = 0.6
farm_options.turbine_options.diameter = 20
farm_options.upwind_correction = True

turbine_location = []
# for i in range(850,1200,200):
#     for j in range(250, 400, 50):
#         turbine_location.append([i,j])
# for i in range(950,1200,200):
#     for j in range(275, 400, 50):
#         turbine_location.append([i,j])

# for i in range(850,1200,100):
#     for j in range(250, 400, 50):
#         turbine_location.append([i,j])

locations = [955.4869483731437, 248.87995382497417, 960.841980078896, 290.65184567303214, 963.1722783852147, 336.9853603396593, 986.5212212062034, 223.6435929578171, 1029.257875136302, 279.73098016183457, 1039.6084355824091, 359.44423976226136, 991.8973452862559, 265.44144520323, 998.202509928964, 304.94138063162, 1002.4176078238502, 344.7186725188866, 1022.9628258262842, 240.13630203637015, 1035.5630397789905, 319.230915589998, 1074.1193249997182, 379.66795138077464]
for i in range(int(len(locations)/2)):
    turbine_location.append([locations[2*i],locations[2*i+1]])

farm_options.turbine_coordinates =[
    [Constant(xy[0]),Constant(xy[1])] for xy in turbine_location
    ]
farm_options.considering_yaw = True

angles = [-15.907106886401639, 8.020416477205213, 18.56072483217048, -32.63483670963606, -23.927349620249906, 35.823709040288215, -20.938312764939862, 16.171681348869846, 28.711942595496982, -44.64497235214166, 25.943879783236305, 38.43039514486168]

farm_options.turbine_axis = [Constant(angle) for angle in angles]

# farm_options.farm_alpha = Function(P1_2d)
#add turbines to SW_equations
options.discrete_tidal_turbine_farms[2] = farm_options

def update_forcings(t):
    print_output("Updating tidal elevation at t = {}".format(t))
    tidal_vel.project(-Constant(speed))

#set initial condition
solver_obj.assign_initial_conditions(uv=as_vector((-Constant(speed), 0.0)), elev=tidal_elev)

# Operation of tidal turbine farm through a callback
cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')

# start computer forward model
solver_obj.iterate(update_forcings=update_forcings)

t_end = time.time()
print('time cost: {0:.2f}min'.format((t_end - t_start)/60))




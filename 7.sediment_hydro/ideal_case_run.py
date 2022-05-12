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


t_start = time.time()

get_index = os.path.basename(sys.argv[0])
the_fac = float(get_index[:-3])

# the_fac = 16
myparameters={}
myparameters['conservative'] = False
myparameters['speed'] = 2 * the_fac/10
myparameters['run_time'] = 1*2
myparameters['morfac'] = 100
myparameters['bath_value'] = 40
myparameters['sediment_size'] = 16*1e-5
myparameters['bnd_sediment'] = 0
myparameters['diffusivity'] = 0.15
myparameters['h_viscosity'] = 1e-6
myparameters['mor_viscosity'] = 1e-6
myparameters['bed_reference_height'] = 3 * myparameters['sediment_size']
myparameters['porosity'] = 0.4
myparameters['bedload'] = True
myparameters['suspended'] = True

# choose directory to output results
st = datetime.datetime.fromtimestamp(t_start).strftime('%Y-%m-%d')# %H:%M:%S')
special_mark = 'inflow_velocity/'+str(float(the_fac))
output_dir = '../../outputs/sediment/results/'+special_mark

### set up the Thetis solver obj as usual ##
mesh2d = Mesh('./mesh.msh')

dt = 30
end_time = 3600*24*myparameters['run_time']

#set viscosity bumps at in-flow boundaries.
P1_2d = FunctionSpace(mesh2d, 'CG', 1)
V = FunctionSpace(mesh2d, "CG", 1)
DG_2d = FunctionSpace(mesh2d, "DG", 1)
vector_dg = VectorFunctionSpace(mesh2d, "DG", 1)


bathymetry_2d = Function(V, name='bathymetry_2d').assign(Constant(myparameters['bath_value']))
# chk = DumbCheckpoint('/home/can/Git_thetis/outputs/sediment/test/'+'M100-H40-24h-2022-03-03-outputs-final_bathymetry',mode=FILE_READ)
# bathymetry_2d = Function(V, name='final_bathymetry')
# chk.load(bathymetry_2d)
# chk.close()

initial_bathymetry_2d = Function(V, name='ini_bathymetry_2d').assign(Constant(myparameters['bath_value']))

# initialise velocity and elevation
chk = DumbCheckpoint("../../outputs/sediment/pure_hydro-v1/elevation", mode=FILE_READ)
elev = Function(DG_2d, name="elevation")
chk.load(elev)
chk.close()

chk = DumbCheckpoint('../../outputs/sediment/pure_hydro-v1/velocity', mode=FILE_READ)
uv = Function(vector_dg, name="velocity")
chk.load(uv)
chk.close()

# create solver and set options
solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry_2d)
options = solver_obj.options

options.sediment_model_options.solve_suspended_sediment = myparameters['suspended']
options.sediment_model_options.use_bedload = myparameters['bedload']
options.sediment_model_options.solve_exner = True
# options.sediment_model_options.use_angle_correction = True
# options.sediment_model_options.use_slope_mag_correction = True
# options.sediment_model_options.use_secondary_current = True
# options.sediment_model_options.use_advective_velocity_correction = False

options.sediment_model_options.use_sediment_conservative_form = myparameters['conservative']
options.sediment_model_options.average_sediment_size = Constant(myparameters['sediment_size'])
options.sediment_model_options.bed_reference_height = Constant(myparameters['bed_reference_height'])
options.sediment_model_options.morphological_acceleration_factor = Constant(myparameters['morfac'])
options.sediment_model_options.morphological_viscosity = Constant(myparameters['mor_viscosity'])
options.sediment_model_options.porosity = Constant(myparameters['porosity'])
options.horizontal_viscosity = Constant(myparameters['h_viscosity'])
options.nikuradse_bed_roughness = Constant(3*options.sediment_model_options.average_sediment_size)
options.horizontal_diffusivity = Constant(myparameters['diffusivity'])

options.simulation_end_time = end_time/myparameters['morfac']
options.simulation_export_time = dt#options.simulation_end_time/45
options.output_directory = output_dir
options.check_volume_conservation_2d = True

if options.sediment_model_options.solve_suspended_sediment:
    options.fields_to_export = ['sediment_2d', 'uv_2d', 'elev_2d', 'bathymetry_2d']  # note exporting bathymetry must be done through export func
    options.sediment_model_options.check_sediment_conservation = True
else:
    options.fields_to_export = ['uv_2d', 'elev_2d', 'bathymetry_2d']  # note exporting bathymetry must be done through export func
options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d', 'bathymetry_2d']
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
options.norm_smoother = Constant(0.1)

# using nikuradse friction


if not hasattr(options.timestepper_options, 'use_automatic_timestep'):
    options.timestep = dt

# assign boundary conditions
left_tag = 1
right_tag = 2
coasts_tag = 3

tidal_elev = Function(get_functionspace(mesh2d, "CG", 1), name='tidal_elev').assign(0.0)


solver_obj.bnd_functions['shallow_water'] = {
    right_tag: {'flux':Constant(-myparameters['speed']*40*300)},
    left_tag: {'elev': tidal_elev},
}

# Initialise Discrete turbine farm characteristics
farm_options = DiscreteTidalTurbineFarmOptions()
farm_options.turbine_type = 'constant'
farm_options.turbine_options.thrust_coefficient = 0.6
farm_options.turbine_options.diameter = 20
farm_options.upwind_correction = True

turbine_location = []

farm_options.turbine_coordinates =[
    [Constant(700),Constant(150)]
    ]
farm_options.considering_yaw = False

#add turbines to SW_equations
options.discrete_tidal_turbine_farms[2] = farm_options

# def update_forcings(t):
#     print_output("Updating tidal elevation at t = {}".format(t))
#     tidal_vel.project(-Constant(speed))

#set initial condition
if options.sediment_model_options.solve_suspended_sediment:
    solver_obj.bnd_functions['sediment'] = {
        right_tag:{'flux':Constant(-myparameters['bnd_sediment']*40*300),'equilibrium':None},
        left_tag:{'elev': tidal_elev}
    }
    solver_obj.assign_initial_conditions(uv=uv, elev=elev)
else:
    solver_obj.assign_initial_conditions(uv=uv, elev=elev)

# Operation of tidal turbine farm through a callback
cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')


# start computer forward model
solver_obj.iterate()#update_forcings=update_forcings)

t_end = time.time()
print('time cost: {0:.2f}min'.format((t_end - t_start)/60))

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank == 0 :
    with open(output_dir+'/myparameters.txt','a+') as f:
        f.write(str(myparameters))
    yag = yagmail.SMTP(user = '623001493@qq.com',password = 'ouehigyjxpidbbcj', host = 'smtp.qq.com')
    yag.send(to = ['623001493@qq.com'], subject = 'Python done', contents = ['College computer Python Finished'])
else:
    pass
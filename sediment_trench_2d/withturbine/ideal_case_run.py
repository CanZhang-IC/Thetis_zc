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
import matplotlib.pyplot as plt

t_start = time.time()

conservative = False

H = 40
distance = 10
speed = 2
output_dir = '../../../outputs/sediment/run'#'../../../outputs/4.yaw/Yaw_Ideal/op-conference_mesh2-5_40/f20-cos00/test'#forward-aligned-both-op'
#Comment for testing forward model

### set up the Thetis solver obj as usual ##
mesh2d = Mesh('../../prepare_ideal_meshes/rectangular.msh')

# timestep = 5
# t_export = 2 * timestep 
# t_end = 6*20*t_export #12000

morfac = 100
dt = 1
end_time = 30*3600

diffusivity = 0.15

#set viscosity bumps at in-flow boundaries.
P1_2d = FunctionSpace(mesh2d, 'CG', 1)
V = FunctionSpace(mesh2d, "CG", 1)
DG_2d = FunctionSpace(mesh2d, "DG", 1)
vector_dg = VectorFunctionSpace(mesh2d, "DG", 1)

h_viscosity = Constant(1e-6)

bathymetry_2d = Function(V, name='bathymetry_2d').assign(Constant(40))

# initialise velocity and elevation
chk = DumbCheckpoint("../../../outputs/sediment/hydrodynamics_trench/elevation", mode=FILE_READ)
elev = Function(DG_2d, name="elevation")
chk.load(elev)
chk.close()

chk = DumbCheckpoint('../../../outputs/sediment/hydrodynamics_trench/velocity', mode=FILE_READ)
uv = Function(vector_dg, name="velocity")
chk.load(uv)
chk.close()

# create solver and set options
solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry_2d)
options = solver_obj.options

options.sediment_model_options.solve_suspended_sediment = True
options.sediment_model_options.use_bedload = True
options.sediment_model_options.solve_exner = True

options.sediment_model_options.use_sediment_conservative_form = conservative
options.sediment_model_options.average_sediment_size = Constant(160*(10**(-6)))
options.sediment_model_options.bed_reference_height = Constant(0.025)
options.sediment_model_options.morphological_acceleration_factor = Constant(morfac)

options.simulation_end_time = end_time/morfac
options.simulation_export_time = options.simulation_end_time/45
options.output_directory = output_dir
options.check_volume_conservation_2d = True

if options.sediment_model_options.solve_suspended_sediment:
    options.fields_to_export = ['sediment_2d', 'uv_2d', 'elev_2d', 'bathymetry_2d']  # note exporting bathymetry must be done through export func
    options.sediment_model_options.check_sediment_conservation = True
else:
    options.fields_to_export = ['uv_2d', 'elev_2d', 'bathymetry_2d']  # note exporting bathymetry must be done through export func

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
options.nikuradse_bed_roughness = Constant(3*options.sediment_model_options.average_sediment_size)
options.horizontal_viscosity = h_viscosity
options.horizontal_diffusivity = Constant(diffusivity)

if not hasattr(options.timestepper_options, 'use_automatic_timestep'):
    options.timestep = dt

# assign boundary conditions
left_tag = 1
right_tag = 2
coasts_tag = 3
tidal_elev = Function(get_functionspace(mesh2d, "CG", 1), name='tidal_elev').assign(0.0)

tidal_vel = Function(V).assign(0.0)

solver_obj.bnd_functions['shallow_water'] = {
    right_tag: {'uv':as_vector((-Constant(speed), 0.0))},
    left_tag: {'elev': tidal_elev},
}

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

for i in range(850,1200,100):
    for j in range(250, 400, 50):
        turbine_location.append([i,j])

farm_options.turbine_coordinates =[
    # [Constant(xy[0]),Constant(xy[1])] for xy in turbine_location
    [Constant(500),Constant(50)]
    ]
farm_options.considering_yaw = False

#add turbines to SW_equations
options.discrete_tidal_turbine_farms[2] = farm_options

# def update_forcings(t):
#     print_output("Updating tidal elevation at t = {}".format(t))
#     tidal_vel.project(-Constant(speed))

#set initial condition
solver_obj.assign_initial_conditions(uv=as_vector((-Constant(speed), 0.0)), elev=tidal_elev)

# # Operation of tidal turbine farm through a callback
# cb = turbines.TurbineFunctionalCallback(solver_obj)
# solver_obj.add_callback(cb, 'timestep')

# start computer forward model
solver_obj.iterate()#update_forcings=update_forcings)

# record final bathymetry for plotting
xaxisthetis1 = []
baththetis1 = []

for i in np.linspace(100, 900, 1000):
    xaxisthetis1.append(i)
    if conservative:
        baththetis1.append(-solver_obj.fields.bathymetry_2d.at([i, 50]))
    else:
        baththetis1.append(-solver_obj.fields.bathymetry_2d.at([i, 50]))

if os.getenv('THETIS_REGRESSION_TEST') is None:
    # Compare model and experimental results
    # (this part is skipped when run as a test)
    # data = np.genfromtxt('experimental_data.csv', delimiter=',')

    # plt.scatter([i[0] for i in data], [i[1] for i in data], label='Experimental Data')

    plt.plot(xaxisthetis1, baththetis1, label='Thetis')
    plt.legend()
    plt.show()




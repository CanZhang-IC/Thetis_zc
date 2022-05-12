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
import rmse_r2
import numpy
import time
import h5py

t_start = time.time()

get_index = os.path.basename(sys.argv[0])
namelength = len('forward')
P_factor = float(get_index[namelength:-3])

# if angle_H[0] == '0':
#     angle = 0
#     H = int(angle_H[1:])
# else:
#     angle = int(angle_H[:2])
#     H = int(angle_H[2:])

angle, H = 0, 40 

speed = 2
output_dir = '../../../outputs/6.yaw_environment/Paper3/Yaw_Ideal-conference_mesh2_with_effected_area/forward/fix-l_y_itc/test/P_factor_'+str(P_factor)+'_op'

### set up the Thetis solver obj as usual ##
mesh2d = Mesh('../../prepare_ideal_meshes/conference_mesh2_with_effected_area.msh')

timestep = 60
t_export = 2 * timestep
t_end = 20*t_export #12000


#set viscosity bumps at in-flow boundaries.
P1_2d = FunctionSpace(mesh2d, 'CG', 1)
x = SpatialCoordinate(mesh2d)
h_viscosity = Constant(1e-3)
# h_viscosity = Function(P1_2d).interpolate(conditional(le(x[0], 20), 20+viscosity-x[0], conditional(ge(x[0],1980),x[0]-viscosity,viscosity)))
# File(output_dir+'/viscosity.pvd').write(h_viscosity)

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
farm_options.turbine_options.C_support = 1
farm_options.turbine_options.A_support = H/2
farm_options.upwind_correction = True

# turbine_location = []

# for i in range(850,1200,200):
#     for j in range(250, 400, 50):
#         turbine_location.append([i,j])
# for i in range(950,1200,200):
#     for j in range(275, 400, 50):
#         turbine_location.append([i,j])

# for i in range(1450,1800,100):
#     for j in range(250, 400, 500):
        # turbine_location.append([i,j])

result_output_dir = '../../../outputs/6.yaw_environment/Paper3/Yaw_Ideal-conference_mesh2_with_effected_area/fix-l_y_itc/P_factor_'+str(P_factor)+'_op'
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank == 0:
    def_file = h5py.File(result_output_dir+'/diagnostic_'+'controls'+'.hdf5','r+')
    for name, data in def_file.items():
        all_controls = list(data[-1])
        iteration_numbers = len(data)
else:
    all_controls = None
    iteration_numbers = None
all_controls = comm.bcast(all_controls,root = 0)
iteration_numbers = comm.bcast(iteration_numbers, root = 0)

all_controls =  [1441.711348941185, 258.8485539883428, 1466.692608441959, 292.803535563495, 1447.838359631565, 343.4952778181821, 1556.7316426836678, 210.0, 1548.7751067425884, 326.3368387278963, 1553.7641700549616, 375.22990450992654, 1648.6960552657824, 227.91759268405687, 1653.8166074610615, 277.31250372432737, 1650.32436484169, 390.0, 1749.7747896872788, 246.14428636712614, 1752.789467864781, 308.92580263731367, 1748.7640637697853, 356.9118445169647] + all_controls

turbine_location = all_controls[:24]
farm_options.turbine_coordinates =[
    # [Constant(xy[0]),Constant(xy[1])] for xy in turbine_location
    [Constant(turbine_location[2*i]),Constant(turbine_location[2*i+1])] for i in range(int(len(turbine_location)/2))
    ]

farm_options.considering_yaw = True
# farm_options.turbine_axis = [Constant(angle) for i in range(len(farm_options.turbine_coordinates))]
angle_outputs = all_controls[24:36]
farm_options.turbine_axis = [Constant(angle) for angle in angle_outputs]

farm_options.considering_individual_thrust_coefficient = True
# farm_options.individual_thrust_coefficient = [Constant(0.6) for i in range(len(farm_options.turbine_coordinates))]
itc_outputs = all_controls[36:]
farm_options.individual_thrust_coefficient = [Constant(itc) for itc in itc_outputs]

# farm_options.farm_alpha = Function(P1_2d)
#add turbines to SW_equations
options.discrete_tidal_turbine_farms[2] = farm_options

def update_forcings(t):
    print_output("Updating tidal elevation at t = {}".format(t))
    tidal_vel.project(-Constant(speed))

solver_obj.assign_initial_conditions(uv=as_vector((-Constant(speed), 0.0)), elev=tidal_elev)

# Operation of tidal turbine farm through a callback
cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')

#Effected area location
E_area_centre_point = [450,250]
E_area_circle = 40

# Operation of tidal turbine farm about each turbine output through a callback
cb2 = rmse_r2.RMSECallback(solver_obj,'../../../outputs/6.yaw_environment/Yaw_Ideal/restart-conference_mesh2_with_effected_area', E_area_centre_point, E_area_circle)
solver_obj.add_callback(cb2,'timestep')

# start computer forward model
solver_obj.iterate(update_forcings=update_forcings)

# ###set up interest functional and control###
power_output= sum(cb.current_power)
maxoutput, maxeffect = 8057.610315283041, 4872.271654375774#10567.41978987577, 6277.669515955061 #
interest_functional = (P_factor*(power_output/maxoutput)-(1-P_factor)*(cb2.RMSE_current[-1]/maxeffect))*maxoutput

if rank ==0:

    with open('result.txt','a+') as f:
        f.write(str(P_factor)+'\t')
        f.write(str(interest_functional)+'\t'+str(power_output)+'\t'+str(cb2.RMSE_current[-1])+'\t')
        f.write(str(iteration_numbers) +'\n')
else:
    pass



t_end = time.time()
print('time cost: {0:.2f}min'.format((t_end - t_start)/60))



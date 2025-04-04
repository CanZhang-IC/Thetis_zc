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
import yagmail

t_start = time.time()

# get_index = os.path.basename(sys.argv[0])
# namelength = len('closed_boundary_upwind')
# angle = get_index[namelength:-3]

angle1,angle2 = 10, 0
H = 0.5
speed = -0.33
output_name = 'mesh_0.05_0.5_f30/angle_'+str(angle1)+'_'+str(angle2)
output_dir = '../../../outputs/4.yaw/Yaw_Ideal/experiments/'+output_name
#Comment for testing forward model

### set up the Thetis solver obj as usual ##
mesh2d = Mesh('../../prepare_ideal_meshes/experiment_mesh.msh')

tidal_amplitude = 5.
tidal_period = 12.42*60*60
timestep = 30
t_export = 2 * timestep
t_end = 20*t_export #12000


#set viscosity bumps at in-flow boundaries.
P1_2d = FunctionSpace(mesh2d, 'CG', 1)
x = SpatialCoordinate(mesh2d)
# h_viscosity = Constant(1e-2)
viscosity = 1e-3
h_viscosity = Function(P1_2d).interpolate(conditional(le(x[0], 1), 1+viscosity-x[0], viscosity))
File(output_dir+'/viscosity.pvd').write(h_viscosity)



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
farm_options.turbine_options.thrust_coefficient = 0.8
farm_options.turbine_options.diameter = 0.3
farm_options.turbine_options.C_support = 1
farm_options.turbine_options.A_support = 0.06*0.5
farm_options.upwind_correction = True

location_turbine = []

# for i in [0,1,2,3,4,6,8,10]:
#     location_turbine.append([25+0.3*i, 2.5])
# farm_options.turbine_coordinates = [[Constant(xy[0]),Constant(xy[1])] for xy in location_turbine]
farm_options.turbine_coordinates =[
    [Constant(25),Constant(2.5)],[Constant(25+5*farm_options.turbine_options.diameter),Constant(2.5)]
    ]

farm_options.considering_yaw = True
farm_options.turbine_axis = [Constant(angle1),Constant(angle2)]

# farm_options.farm_alpha = Function(P1_2d)
#add turbines to SW_equations
options.discrete_tidal_turbine_farms[2] = farm_options

def update_forcings(t):
    print_output("Updating tidal elevation at t = {}".format(t))
    tidal_vel.project(-Constant(speed))

    uv, eta = split(solver_obj.fields.solution_2d)
    # alpha_flood = sum((conditional((x[0]-xi)**2+(x[1]-yi)**2 < farm_options.turbine_options.diameter**2, alphai/180*pi, 0) \
    #             for alphai,(xi,yi) in zip(farm_options.turbine_axis[:len(farm_options.turbine_coordinates)],farm_options.turbine_coordinates)))
    # alpha_ebb   = sum((conditional((x[0]-xi)**2+(x[1]-yi)**2 < farm_options.turbine_options.diameter**2, alphai/180*pi, 0) \
    #             for alphai,(xi,yi) in zip(farm_options.turbine_axis[len(farm_options.turbine_coordinates):],farm_options.turbine_coordinates)))
    # alphainterpolation = conditional(uv[0] > 0, alpha_flood, alpha_ebb)
    # farm_options.farm_alpha.interpolate(alphainterpolation)

#set initial condition
solver_obj.assign_initial_conditions(uv=as_vector((-Constant(speed), 0.0)), elev=tidal_elev)

# Operation of tidal turbine farm through a callback
cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')

# Operation of tidal turbine farm about each turbine output through a callback
# cb2 = turbines.EachTurbineFunctionalCallback(solver_obj)
# solver_obj.add_callback(cb2, 'timestep')

#detector location in the wake area
locations=[]
names=[]
xx=[1,2,3,4,6,8,10,12]
yy=[-1,-0.8,-0.6,-0.5,-0.4,-0.2,-0.1, 0, 0.1,0.2,0.4,0.5,0.6,0.8,1]
for i in xx:
    for j in yy:
        locations.append((25+0.3*i,2.5+0.3*j))
        names.append('name'+str((i,j)))
        
cb2 = DetectorsCallback(solver_obj, locations, ['elev_2d', 'uv_2d'], name='detectors',detector_names=names)
solver_obj.add_callback(cb2, 'timestep')

# start computer forward model
solver_obj.iterate(update_forcings=update_forcings)

t_end = time.time()
print('time cost: {0:.2f}min'.format((t_end - t_start)/60))

# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# if rank == 0 :
#     yag = yagmail.SMTP(user = '623001493@qq.com',password = 'Zc623oo1493', host = 'smtp.qq.com')
#     yag.send(to = ['canzhang2019@gmail.com'], subject = 'Python done', contents = [output_name])
# else:
#     pass

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
import h5py

t_start = time.time()

H = 40
output_dir = '../../../outputs/4.yaw/Conference2/test_v0.1_forward_angle'
#Comment for testing forward model
test_gradient = False
optimise = False

### set up the Thetis solver obj as usual ###
mesh2d = Mesh('../../prepare_ideal_meshes/headland2.msh')

tidal_amplitude = 5.
tidal_period = 12.42*60*60
timestep = 200
t_export = 600
t_end = 100*t_export + tidal_period


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
# options.timestepper_options.use_semi_implicit_linearization = True
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

# Initialise Discrete turbine farm characteristics
farm_options = DiscreteTidalTurbineFarmOptions()
farm_options.turbine_type = 'constant'
farm_options.turbine_options.thrust_coefficient = 0.8
farm_options.turbine_options.diameter = 20
farm_options.upwind_correction = False

site_x1, site_x2, site_y1, site_y2 = 920, 1080, 250, 350
turbinelocation = []
for i in range(940,1080,40):
    for j in range(275,350,26):
        turbinelocation.append([i,j])
farm_options.turbine_coordinates =[[Constant(i[0]),Constant(i[1])] for i in turbinelocation]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0 :
    # location = []
    angle = []
    df = h5py.File('../../../outputs/4.yaw/Conference2/yaw_aligned/diagnostic_controls.hdf5','r+')
    for name,data in df.items():
        for i in range(24):
            # location.append(data[-1][i])
            angle.append(data[-1][i])
    print(angle)   
else:
    # location = None
    angle = None
# location = comm.bcast(location, root=0)
angle = comm.bcast(angle,root = 0)

# farm_options.turbine_coordinates =[[Constant(location[2*i]),Constant(location[2*i+1])] for i in range(12)]

farm_options.considering_yaw = True
farm_options.turbine_axis = [Constant(i) for i in angle]
farm_options.farm_alpha = Function(P1_2d)
#add turbines to SW_equations
options.discrete_tidal_turbine_farms[2] = farm_options


def update_forcings(t):
    print_output("Updating tidal elevation at t = {}".format(t))
    tidal_elev.project(tidal_amplitude*sin(omega*t + omega/Constant(pow(g*H, 0.5))*x[0]))

    uv, eta = split(solver_obj.fields.solution_2d)
    alpha_flood = sum((conditional((x[0]-xi)**2+(x[1]-yi)**2 < farm_options.turbine_options.diameter**2, alphai/180*pi, 0) \
                for alphai,(xi,yi) in zip(farm_options.turbine_axis[:len(farm_options.turbine_coordinates)],farm_options.turbine_coordinates)))
    alpha_ebb = sum((conditional((x[0]-xi)**2+(x[1]-yi)**2 < farm_options.turbine_options.diameter**2, alphai/180*pi, 0) \
        for alphai,(xi,yi) in zip(farm_options.turbine_axis[len(farm_options.turbine_coordinates):],farm_options.turbine_coordinates)))
    alphainterpolation = conditional(uv[0] > 0, alpha_flood, alpha_ebb)
    farm_options.farm_alpha.interpolate(alphainterpolation)

#set initial condition
# solver_obj.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=tidal_elev)
solver_obj.load_state(100, '../../../outputs/4.yaw/Conference2/restart')

# Operation of tidal turbine farm through a callback
cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')

cb2 = turbines.EachTurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb2,'timestep')

# start computer forward model
solver_obj.iterate(update_forcings=update_forcings)

t_end = time.time()
print('time cost: {0:.2f}min'.format((t_end - t_start)/60))

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0 :
    yag = yagmail.SMTP(user = '623001493@qq.com',password = 'lovezhongqiu', host = 'smtp.qq.com')
    yag.send(to = ['canzhang2019@gmail.com'], subject = 'Python done', contents = ['Forward run finishes. The time cost is:' + str(int((t_end-t_start)/60)) + 'min.'])
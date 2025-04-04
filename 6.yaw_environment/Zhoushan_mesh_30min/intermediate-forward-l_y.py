from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
import numpy
op2.init(log_level=INFO)
import sys
sys.path.append('../..')
import prepare.utm, prepare.myboundary_30min, prepare.detectors
import os
import time
import yagmail
import rmse_r2
import h5py

t_start = time.time()

file_dir = '../../'

# get_index = os.path.basename(sys.argv[0])
# namelength = len('intermediate-forward-l_y')
# P_factor = float(get_index[namelength:-3])

P_factor = 1

output_dir = '../../../outputs/6.yaw_environment/Paper3/Zhoushan_mesh_30min/optimisation/intermediate-forward-l_y-P_factor_'+str(P_factor)


mesh2d = Mesh(file_dir+'mesh/mesh.msh')

#timestepping options
dt = 30*60 # reduce this if solver does not converge
t_export = 30*60 
# t_end = 1555200
# t_end = 1216800+ 13*60*60 # spring
# t_end = 885600 + 13*60*60 # middle
# t_end = 612000 + 13*60*60 # neap
t_end = 898200 + 5*60*60 # middle

P1 = FunctionSpace(mesh2d, "CG", 1)

# read bathymetry code
chk = DumbCheckpoint(file_dir+'prepare/bathymetry', mode=FILE_READ)
bathymetry2d = Function(P1)
chk.load(bathymetry2d, name='bathymetry')
chk.close()

#read viscosity / manning boundaries code
chk = DumbCheckpoint(file_dir+'prepare/viscosity', mode=FILE_READ)
h_viscosity = Function(P1, name='viscosity') 
chk.load(h_viscosity)
chk.close()

#manning = Function(P1,name='manning')
chk = DumbCheckpoint(file_dir+'prepare/manning', mode=FILE_READ)
manning = Function(bathymetry2d.function_space(), name='manning')
chk.load(manning)
chk.close()

#account for Coriolis code
def coriolis(mesh, lat,):
    R = 6371e3
    Omega = 7.292e-5
    lat_r = lat * pi / 180.
    f0 = 2 * Omega * sin(lat_r)
    beta = (1 / R) * 2 * Omega * cos(lat_r)
    x = SpatialCoordinate(mesh)
    x_0, y_0, utm_zone, zone_letter = prepare.utm.from_latlon(lat, 0)
    coriolis_2d = Function(FunctionSpace(mesh, 'CG', 1), name="coriolis_2d")
    coriolis_2d.interpolate(f0 + beta * (x[1] - y_0))
    return coriolis_2d
coriolis_2d =coriolis(mesh2d, 30)



# --- create solver ---
solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry2d)
options = solver_obj.options
options.cfl_2d = 1.0
options.use_nonlinear_equations = True
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.coriolis_frequency = coriolis_2d
options.output_directory = output_dir
options.check_volume_conservation_2d = True
options.fields_to_export = ['uv_2d', 'elev_2d']
options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']
options.element_family = "dg-dg"
options.timestepper_type = 'CrankNicolson'
options.timestepper_options.implicitness_theta = 1 #Implicitness parameter, default 0.5
options.timestepper_options.use_semi_implicit_linearization = True#If True use a linearized semi-implicit scheme
options.use_wetting_and_drying = True #Wetting and drying is included through the modified bathymetry formulation of Karna et al. (2011). 
options.wetting_and_drying_alpha = Constant(0.5) #need to check if this is a good value
options.manning_drag_coefficient = manning 
#options.quadratic_drag_coefficient = Constant(0.0025)
options.horizontal_viscosity = h_viscosity #the viscosity 'cushion' we created in initialisation & loaded above
#Epshteyn and Riviere (2007). Estimation of penalty parameters for symmetric interior penalty Galerkin methods. Journal of Computational and Applied Mathematics, 206(2):843-872. http://dx.doi.org/10.1016/j.cam.2006.08.029
options.use_grad_div_viscosity_term = True
options.use_grad_depth_viscosity_term = False
options.timestep = dt
options.timestepper_options.solver_parameters = {'snes_monitor': None,
                                                 'snes_rtol': 1e-5,
                                                 'ksp_type': 'preonly',
                                                 'pc_type': 'lu',
                                                 'pc_factor_mat_solver_type': 'mumps',
                                                 'mat_type': 'aij'
                                                 }

# set boundary/initial conditions code
tidal_elev = Function(bathymetry2d.function_space())
tidal_v = Function(VectorFunctionSpace(mesh2d,"CG",1))
solver_obj.bnd_functions['shallow_water'] = {
        200: {'elev': tidal_elev,'uv': tidal_v},  #set open boundaries to tidal_elev function
  }


# Initialise Discrete turbine farm characteristics
farm_options = DiscreteTidalTurbineFarmOptions()
farm_options.turbine_type = 'constant'
farm_options.turbine_options.thrust_coefficient = 0.6
farm_options.turbine_options.diameter = 20
# farm_options.turbine_options.C_support = 1
# farm_options.turbine_options.A_support = H/2
farm_options.upwind_correction = True

xmin,ymin,xmax,ymax = 443340, 3322634, 443592, 3322848 

result_output_dir = '../../../outputs/6.yaw_environment/Paper3/Zhoushan_mesh_30min/optimisation/intermediate-op-l_y_vary-P_factor_'+str(P_factor)
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


turbine_location = all_controls[:24]

farm_options.turbine_coordinates =[[Constant(turbine_location[2*i]),Constant(turbine_location[2*i+1])] for i in range(int(len(turbine_location)/2))]

farm_options.considering_yaw = True
flood_dir,ebb_dir = all_controls[24:36], all_controls[36:]

farm_options.turbine_axis = [Constant(i) for i in flood_dir] + [Constant(i) for i in ebb_dir]

farm_options.considering_individual_thrust_coefficient = False
farm_options.individual_thrust_coefficient = [Constant(0.6) for i in range(len(farm_options.turbine_axis))]

#add turbines to SW_equations
options.discrete_tidal_turbine_farms[2] = farm_options

def update_forcings(t):
  with timed_stage('update forcings'):
    print_output("Updating tidal field at t={}".format(t))
    elev = prepare.myboundary_30min.set_tidal_field(Function(bathymetry2d.function_space()), t, dt)
    tidal_elev.project(elev) 
    v = prepare.myboundary_30min.set_velocity_field(Function(VectorFunctionSpace(mesh2d,"CG",1)),t,dt)
    tidal_v.project(v)
    print_output("Done updating tidal field")

###spring:676,middle:492,neap:340###
# solver_obj.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=Constant(0.0))
solver_obj.load_state(499, outputdir='../../../outputs/6.yaw_environment/Paper3/Zhoushan_mesh_30min/restart_30min-e&v')

# Operation of tidal turbine farm through a callback
cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')

#Effected area location
E_area_centre_point = [(xmin+xmax)/2,ymax+800]
E_area_circle = 60

# Operation of tidal turbine farm about each turbine output through a callback
cb2 = rmse_r2.RMSECallback(solver_obj,'../../../outputs/6.yaw_environment/Paper3/Zhoushan_mesh_30min/optimisation/intermediate-forward', E_area_centre_point, E_area_circle)
solver_obj.add_callback(cb2,'timestep')

# start computer forward model
solver_obj.iterate(update_forcings=update_forcings)

# ###set up interest functional and control###
power_output= sum(cb.average_power)
maxoutput, maxeffect = 2733.6082917429912, 1877.8060067867286
interest_functional = (P_factor*(power_output/maxoutput)-(1-P_factor)*(cb2.RMSEaverage/maxeffect))*maxoutput
print(interest_functional,power_output,cb2.RMSEaverage)

if rank ==0:

    with open('result.txt','a+') as f:
        f.write(str(P_factor)+'\t')
        f.write(str(interest_functional)+'\t'+str(power_output)+'\t'+str(cb2.RMSEaverage)+'\t')
        f.write(str(iteration_numbers) +'\n')
else:
    pass


t_end = time.time()
print('time cost: {0:.2f}min'.format((t_end - t_start)/60))

if P_factor == 0.9:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0 :
        yag = yagmail.SMTP(user = '623001493@qq.com',password = 'Zc623oo1493', host = 'smtp.qq.com')
        yag.send(to = ['canzhang2019@gmail.com'], subject = 'Python done', contents = ['Location and Yaw Optimisation Finished.'])
    else:
        pass

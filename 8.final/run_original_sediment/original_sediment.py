from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
op2.init(log_level=INFO)
import sys
sys.path.append('../..')
import prepare_vs_otf_oe1.utm, prepare_vs_otf_oe1.myboundary
import prepare_cable.Hybrid_Code
from prepare_cable.cable_overloaded import cablelength
import numpy
import time
import yagmail
import h5py
import rmse_sediment

start_time = time.time()

# get_index = os.path.basename(sys.argv[0])
# BE = float(get_index[:-3])


output_dir = '../../../outputs/8.final/486-original_sediment'
print_output(output_dir[17:])

file_dir = '../../'
mesh2d = Mesh(file_dir+'mesh-vs_otf_oe1/mesh.msh')

###spring:676,middle:492,neap:340###
start_time_point = 486
#timestepping options
dt = 5*60 # reduce this if solver does not converge
t_export = 30*60 
# t_end = 1555200
t_end = 492*t_export+ 301*dt

P1 = FunctionSpace(mesh2d, "CG", 1)

# read bathymetry code
chk = DumbCheckpoint(file_dir+'prepare_vs_otf_oe1/bathymetry', mode=FILE_READ)
bathymetry2d = Function(P1)
chk.load(bathymetry2d, name='bathymetry')
chk.close()

#read viscosity / manning boundaries code
chk = DumbCheckpoint(file_dir+'prepare_vs_otf_oe1/viscosity', mode=FILE_READ)
h_viscosity = Function(P1, name='viscosity')
chk.load(h_viscosity)
chk.close()

#manning = Function(P1,name='manning')
chk = DumbCheckpoint(file_dir+'prepare_vs_otf_oe1/manning', mode=FILE_READ)
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
    x_0, y_0, utm_zone, zone_letter = prepare_vs_otf_oe1.utm.from_latlon(lat, 0)
    coriolis_2d = Function(FunctionSpace(mesh, 'CG', 1), name="coriolis_2d")
    coriolis_2d.interpolate(f0 + beta * (x[1] - y_0))
    return coriolis_2d
coriolis_2d =coriolis(mesh2d, 30)

# --- create solver ---
solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry2d)
options = solver_obj.options

sediment_size,morfac,diffusivity,morphological_viscosity = 40*1e-6,1,0.01,1e-6
options.sediment_model_options.c_bstar_constant = 15/1000
options.sediment_model_options.solve_suspended_sediment = True
options.sediment_model_options.use_bedload = False
options.sediment_model_options.solve_exner = False
# options.sediment_model_options.use_angle_correction = True
# options.sediment_model_options.use_slope_mag_correction = True
# options.sediment_model_options.use_secondary_current = True
# options.sediment_model_options.use_advective_velocity_correction = False

options.sediment_model_options.use_sediment_conservative_form = True
options.sediment_model_options.average_sediment_size = sediment_size
options.sediment_model_options.bed_reference_height = 4*sediment_size
options.sediment_model_options.morphological_acceleration_factor = Constant(morfac)
options.sediment_model_options.morphological_viscosity = Constant(morphological_viscosity)
options.sediment_model_options.porosity = Constant(0.4)
options.horizontal_viscosity = h_viscosity
# define parameters
options.nikuradse_bed_roughness = Constant(3*sediment_size)
options.horizontal_diffusivity = Constant(diffusivity)
options.norm_smoother = Constant(1)

options.cfl_2d = 1.0
options.use_nonlinear_equations = True
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.coriolis_frequency = coriolis_2d
options.output_directory = output_dir
options.check_volume_conservation_2d = True
if options.sediment_model_options.solve_suspended_sediment:
    options.fields_to_export = ['sediment_2d', 'uv_2d', 'elev_2d']  # note exporting bathymetry must be done through export func
    options.fields_to_export_hdf5 = ['sediment_2d', 'uv_2d', 'elev_2d']
    options.sediment_model_options.check_sediment_conservation = True
else:
    options.fields_to_export = ['uv_2d', 'elev_2d']  # note exporting bathymetry must be done through export func
    options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']
options.element_family = "dg-dg"
options.timestepper_type = 'CrankNicolson'
options.timestepper_options.implicitness_theta = 1 #Implicitness parameter, default 0.5
options.timestepper_options.use_semi_implicit_linearization = True#If True use a linearized semi-implicit scheme
options.use_wetting_and_drying = True #Wetting and drying is included through the modified bathymetry formulation of Karna et al. (2011). 
options.wetting_and_drying_alpha = Constant(0.5) #need to check if this is a good value
# options.manning_drag_coefficient = manning 
#options.quadratic_drag_coefficient = Constant(0.0025)
# options.horizontal_viscosity = h_viscosity #the viscosity 'cushion' we created in initialisation & loaded above
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
        200: {'elev': tidal_elev},  #set open boundaries to tidal_elev function
  }


# Initialise Discrete turbine farm characteristics
farm_options = DiscreteTidalTurbineFarmOptions()
farm_options.turbine_type = 'constant'
farm_options.turbine_options.thrust_coefficient = 0.6
farm_options.turbine_options.diameter = 20
farm_options.upwind_correction = True

xmin,ymin,xmax,ymax = 443100, 3322600, 443600, 3323100 
d = farm_options.turbine_options.diameter
farm_options.considering_yaw =  False

result_output_dir = '../../../outputs/8.final/just_power'
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

farm_options.turbine_coordinates = [[Constant(all_controls[2*i]),Constant(all_controls[2*i+1])] for i in range(64)]

#add turbines to SW_equations
options.discrete_tidal_turbine_farms[2] = farm_options

def update_forcings(t):
  with timed_stage('update forcings'):
    print_output("Updating tidal field at t={}".format(t))
    elev = prepare_vs_otf_oe1.myboundary.set_tidal_field(Function(bathymetry2d.function_space()), t, dt)
    tidal_elev.project(elev) 
    v = prepare_vs_otf_oe1.myboundary.set_velocity_field(Function(VectorFunctionSpace(mesh2d,"CG",1)),t,dt)
    tidal_v.project(v)
    print_output("Done updating tidal field")

###spring:676,middle:492,neap:340###
solver_obj.load_state(start_time_point, outputdir='../../../outputs/0.validation/vsotf-discrete-4cores')
#solver_obj.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=Constant(0.0))

# Operation of tidal turbine farm through a callback
cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')

### 考虑环境影响 ###
#Effected area location
E_area_centre_point = [[444500,3319900]] #behind
E_area_circle = [200]
# Operation of tidal turbine farm about each turbine output through a callback
cb2 = rmse_sediment.RMSECallback(solver_obj,'../../../outputs/3.environment/zhoushan-continuous/486-d5-4', E_area_centre_point, E_area_circle,492)
solver_obj.add_callback(cb2,eval_interval='export')

#Effected area location
E_area_centre_point = [[442700,3323760]] #inside
E_area_circle = [100]

# Operation of tidal turbine farm about each turbine output through a callback
cb3 = rmse_sediment.RMSECallback(solver_obj,'../../../outputs/3.environment/zhoushan-continuous/486-d5-4', E_area_centre_point, E_area_circle,492)
solver_obj.add_callback(cb3,eval_interval='export')

#Effected area location
E_area_centre_point1 = [[442973,3321135]] #next to
E_area_circle1 = [200]

# Operation of tidal turbine farm about each turbine output through a callback
cb4 = rmse_sediment.RMSECallback(solver_obj,'../../../outputs/3.environment/zhoushan-continuous/486-d5-4', E_area_centre_point1, E_area_circle1,492)
solver_obj.add_callback(cb4,eval_interval='export')

# start computer forward model
solver_obj.iterate(update_forcings=update_forcings)

###电缆长度计算###
turbine_locations = [float(x) for xy in farm_options.turbine_coordinates for x in xy]
landpointlocation = [444000,3323000]
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank == 0:
    cableclass = prepare_cable.Hybrid_Code.CableCostGA(turbine_locations=turbine_locations, substation_location=landpointlocation,capacity = 6)
    order_w = cableclass.compute_cable_cost_order()
else:
    order_w = []
order_w = comm.bcast(order_w, root=0)



###set up interest functional and control###
### 环境影响
effect_two = sum([cb2.RMSEaverage,cb3.RMSEaverage,cb4.RMSEaverage])
### 电缆长度
landpointlocation_con = [Constant(x) for x in landpointlocation]
order_con = [Constant(i) for j in order_w for i in j]
cablecost = cablelength([x for xy in farm_options.turbine_coordinates for x in xy],landpointlocation_con,order_con)
### 能量输出
power_output= sum(cb.average_power)
### 目标函数组成
interest_functional = power_output

### 水轮机位置信息出图
turbine_density = Function(solver_obj.function_spaces.P1_2d, name='turbine_density')
turbine_density.interpolate(solver_obj.tidal_farms[0].turbine_density)
File(output_dir+'/turbinedensity.pvd').write(turbine_density)

end_time = time.time()
print_output('time cost: {0:.2f}h'.format((end_time - start_time)/60/60))

if 1:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0 :
        yag = yagmail.SMTP(user = '623001493@qq.com',password = 'ouehigyjxpidbbcj', host = 'smtp.qq.com')
        yag.send(to = ['623001493@qq.com'], subject = output_dir[17:], contents = ['Time cose: {0:.2f}h.'.format((end_time-start_time)/60/60)])
    else:
        pass
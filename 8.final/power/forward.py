from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
import numpy
op2.init(log_level=INFO)
import sys
sys.path.append('../..')
import prepare.utm, prepare.myboundary, prepare.detectors
import os
import time
import yagmail
import rmse_r2
import h5py

start_time = time.time()

file_dir = '../../'

# get_index = os.path.basename(sys.argv[0])
# break_point = get_index.index('-')
# xspacing = float(get_index[:break_point])
# yspacing = float(get_index[break_point+1:-3])
# print(str(xspacing)[:-2]+'_'+str(yspacing)[:-2])

tide_period = 'middle'
if tide_period == 'neap':
    t1_end, t2_end = 353, 366
elif tide_period == 'middle':
    t1_end, t2_end = 502, 514
else:
    t1_end, t2_end = 687, 700
xspacing,yspacing = 40.0,40.0

output_dir = '../../../outputs/8.final/flood_ebb-forward/'+str(tide_period)+'-xspacing'+str(xspacing)[:-2] + '-yspacing'+ str(yspacing)[:-2] 
print_output(output_dir[17:])
mesh2d = Mesh(file_dir+'mesh/mesh.msh')

#timestepping options
dt = 5*60 # reduce this if solver does not converge
t_export = 30*60 

t_end1 = 612000 + 13*60*60 # neap
t_end2 = 885600 + 13*60*60 # middle
t_end3 = 1216800+ 13*60*60 # spring

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
options.simulation_end_time = t_end1
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

# --- create solver ---
solver_obj2 = solver2d.FlowSolver2d(mesh2d, bathymetry2d)
options2 = solver_obj2.options
options2.cfl_2d = 1.0
options2.use_nonlinear_equations = True
options2.simulation_export_time = t_export
options2.simulation_end_time = t_end2
options2.coriolis_frequency = coriolis_2d
options2.output_directory = output_dir
options2.check_volume_conservation_2d = True
options2.fields_to_export = ['uv_2d', 'elev_2d']
options2.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']
options2.element_family = "dg-dg"
options2.timestepper_type = 'CrankNicolson'
options2.timestepper_options.implicitness_theta = 1 #Implicitness parameter, default 0.5
options2.timestepper_options.use_semi_implicit_linearization = True#If True use a linearized semi-implicit scheme
options2.use_wetting_and_drying = True #Wetting and drying is included through the modified bathymetry formulation of Karna et al. (2011). 
options2.wetting_and_drying_alpha = Constant(0.5) #need to check if this is a good value
options2.manning_drag_coefficient = manning 
#options2.quadratic_drag_coefficient = Constant(0.0025)
options2.horizontal_viscosity = h_viscosity #the viscosity 'cushion' we created in initialisation & loaded above
#Epshteyn and Riviere (2007). Estimation of penalty parameters for symmetric interior penalty Galerkin methods. Journal of Computational and Applied Mathematics, 206(2):843-872. http://dx.doi.org/10.1016/j.cam.2006.08.029
options2.use_grad_div_viscosity_term = True
options2.use_grad_depth_viscosity_term = False
options2.timestep = dt
options2.timestepper_options.solver_parameters = {'snes_monitor': None,
                                                 'snes_rtol': 1e-5,
                                                 'ksp_type': 'preonly',
                                                 'pc_type': 'lu',
                                                 'pc_factor_mat_solver_type': 'mumps',
                                                 'mat_type': 'aij'
                                                 }

# set boundary/initial conditions code
solver_obj2.bnd_functions['shallow_water'] = {
        200: {'elev': tidal_elev,'uv': tidal_v},  #set open boundaries to tidal_elev function
  }

# --- create solver ---
solver_obj3 = solver2d.FlowSolver2d(mesh2d, bathymetry2d)
options3 = solver_obj3.options
options3.cfl_2d = 1.0
options3.use_nonlinear_equations = True
options3.simulation_export_time = t_export
options3.simulation_end_time = t_end3
options3.coriolis_frequency = coriolis_2d
options3.output_directory = output_dir
options3.check_volume_conservation_2d = True
options3.fields_to_export = ['uv_2d', 'elev_2d']
options3.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']
options3.element_family = "dg-dg"
options3.timestepper_type = 'CrankNicolson'
options3.timestepper_options.implicitness_theta = 1 #Implicitness parameter, default 0.5
options3.timestepper_options.use_semi_implicit_linearization = True#If True use a linearized semi-implicit scheme
options3.use_wetting_and_drying = True #Wetting and drying is included through the modified bathymetry formulation of Karna et al. (2011). 
options3.wetting_and_drying_alpha = Constant(0.5) #need to check if this is a good value
options3.manning_drag_coefficient = manning 
#options3.quadratic_drag_coefficient = Constant(0.0025)
options3.horizontal_viscosity = h_viscosity #the viscosity 'cushion' we created in initialisation & loaded above
#Epshteyn and Riviere (2007). Estimation of penalty parameters for symmetric interior penalty Galerkin methods. Journal of Computational and Applied Mathematics, 206(2):843-872. http://dx.doi.org/10.1016/j.cam.2006.08.029
options3.use_grad_div_viscosity_term = True
options3.use_grad_depth_viscosity_term = False
options3.timestep = dt
options3.timestepper_options.solver_parameters = {'snes_monitor': None,
                                                 'snes_rtol': 1e-5,
                                                 'ksp_type': 'preonly',
                                                 'pc_type': 'lu',
                                                 'pc_factor_mat_solver_type': 'mumps',
                                                 'mat_type': 'aij'
                                                 }

# set boundary/initial conditions code
solver_obj3.bnd_functions['shallow_water'] = {
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

xmin,ymin,xmax,ymax = 443100, 3322600, 443600, 3323100 

# turbine_location = [[443285.232752254, 3322701.82872359], [443251.030737921, 3322795.79798567], [443216.828723589, 3322889.76724775], [443322.820457085, 3322715.50952932], [443288.618442753, 3322809.4787914], [443254.41642842, 3322903.44805348], [443360.408161917, 3322729.19033505], [443326.206147584, 3322823.15959713], [443292.004133252, 3322917.12885921], [443397.995866748, 3322742.87114079], [443363.793852416, 3322836.84040287], [443329.591838083, 3322930.80966495], [443435.58357158, 3322756.55194652], [443401.381557247, 3322850.5212086], [443367.179542915, 3322944.49047068], [443473.171276411, 3322770.23275225], [443438.969262079, 3322864.20201433], [443404.767247746, 3322958.17127641]]
# farm_options.turbine_coordinates =[[Constant(xy[0]),Constant(xy[1])] for xy in turbine_location]

# farm_options.considering_yaw = True


# farm_options.turbine_axis = [Constant(90) for i in range(len(farm_options.turbine_coordinates))] + [Constant(270) for i in range(len(farm_options.turbine_coordinates))]

result_output_dir = '../../../outputs/8.final/flood_ebb/'+str(tide_period)+'-xspacing'+str(xspacing)[:-2] + '-yspacing'+ str(yspacing)[:-2] 
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

farm_options.turbine_coordinates = [[Constant(all_controls[2*i]),Constant(all_controls[2*i+1])] for i in range(18)]

farm_options.considering_yaw = True
flood_dir,ebb_dir = all_controls[36:54], all_controls[54:]

farm_options.turbine_axis = [Constant(i) for i in flood_dir] + [Constant(i) for i in ebb_dir]


#add turbines to SW_equations
options.discrete_tidal_turbine_farms[2] = farm_options
options2.discrete_tidal_turbine_farms[2] = options.discrete_tidal_turbine_farms[2] 
options3.discrete_tidal_turbine_farms[2] = options.discrete_tidal_turbine_farms[2] 

def update_forcings(t):
  with timed_stage('update forcings'):
    print_output("Updating tidal field at t={}".format(t))
    elev = prepare.myboundary.set_tidal_field(Function(bathymetry2d.function_space()), t, dt)
    tidal_elev.project(elev) 
    v = prepare.myboundary.set_velocity_field(Function(VectorFunctionSpace(mesh2d,"CG",1)),t,dt)
    tidal_v.project(v)
    print_output("Done updating tidal field")

###spring:676,middle:492,neap:340###
# solver_obj.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=Constant(0.0))
solver_obj.load_state(340, outputdir='../../../outputs/0.validation/discrete-4cores')

# Operation of tidal turbine farm through a callback
cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')

# start computer forward model
solver_obj.iterate(update_forcings=update_forcings)

###set up interest functional and control###
power_output= sum(cb.average_power)

###spring:676,middle:492,neap:340###
# solver_obj2.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=Constant(0.0))
solver_obj2.load_state(492, outputdir='../../../outputs/0.validation/discrete-4cores')

# Operation of tidal turbine farm through a callback
cb2 = turbines.TurbineFunctionalCallback(solver_obj2)
solver_obj2.add_callback(cb2, 'timestep')

# start computer forward model
solver_obj2.iterate(update_forcings=update_forcings)

###set up interest functional and control###
power_output2= sum(cb2.average_power)

###spring:676,middle:492,neap:340###
# solver_obj3.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=Constant(0.0))
solver_obj3.load_state(676, outputdir='../../../outputs/0.validation/discrete-4cores')

# Operation of tidal turbine farm through a callback
cb3 = turbines.TurbineFunctionalCallback(solver_obj3)
solver_obj3.add_callback(cb3, 'timestep')

# start computer forward model
solver_obj3.iterate(update_forcings=update_forcings)

###set up interest functional and control###
power_output3= sum(cb3.average_power)




power_outputs = [power_output,power_output2,power_output3]
interest_functional = sum(power_outputs)/3

if rank ==0:
    with open('result-l_y.txt','a+') as f:
        f.write(tide_period + '\t' + str(xspacing) + '-' + str(yspacing) + '\t' )
        f.write(str(interest_functional)+'\t'+str(power_outputs)+'\t')
        f.write(str(iteration_numbers) +'\n')
else:
    pass



end_time = time.time()
print_output('time cost: {0:.2f}h'.format((end_time - start_time)/60/60))

if 0:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0 :
        yag = yagmail.SMTP(user = '623001493@qq.com',password = 'ouehigyjxpidbbcj', host = 'smtp.qq.com')
        yag.send(to = ['623001493@qq.com'], subject = output_dir[17:], contents = ['Time cose: {0:.2f}h.'.format((end_time-start_time)/60/60)])
    else:
        pass

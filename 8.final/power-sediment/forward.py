from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
import h5py
op2.init(log_level=INFO)
import sys
sys.path.append('../..')
import prepare.utm, prepare.myboundary, prepare.detectors
import os
import time
import yagmail
import rmse_r2,rmse_r3
import prepare_cable.Hybrid_Code
from prepare_cable.cable_overloaded import cablelength

start_time = time.time()

# get_index = os.path.basename(sys.argv[0])
# break_point = get_index.index('-')
# BEcable = float(get_index[:break_point])
# BEsediment = float(get_index[break_point+1:-3])
# print(str(BEcable)[:-2]+'_'+str(BEsediment)[:-2])
BEcable,BEsediment = 0.0,0.0

###spring:676,middle:492,neap:340###
start_time_point = 492

output_dir = '../../../outputs/8.final/cable-sediment-forward-uv/'+str(BEcable)[:-2] + '-'+ str(BEsediment)[:-2] 
print_output(output_dir[17:])

file_dir = '../../'
mesh2d = Mesh(file_dir+'mesh/mesh.msh')
#timestepping options
dt = 5*60 # reduce this if solver does not converge
t_export = 30*60 
# t_end = 1555200
t_end = 492*t_export+ 13*60*60

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
farm_options.upwind_correction = True


result_output_dir = '../../../outputs/8.final/cable-sediment/'+str(BEcable)[:-2] + '-'+ str(BEsediment)[:-2] 

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

farm_options.turbine_coordinates = [[Constant(all_controls[2*i]),Constant(all_controls[2*i+1])] for i in range(12)]

farm_options.considering_yaw = True
flood_dir,ebb_dir = all_controls[24:36], all_controls[36:]

farm_options.turbine_axis = [Constant(i) for i in flood_dir] + [Constant(i) for i in ebb_dir]

farm_options.considering_individual_thrust_coefficient = False
#add turbines to SW_equations
options.discrete_tidal_turbine_farms[2] = farm_options


def update_forcings(t):
  with timed_stage('update forcings'):
    print_output("Updating tidal field at t={}".format(t))
    elev = prepare.myboundary.set_tidal_field(Function(bathymetry2d.function_space()), t, dt)
    tidal_elev.project(elev) 
    v = prepare.myboundary.set_velocity_field(Function(VectorFunctionSpace(mesh2d,"CG",1)),t,dt)
    tidal_v.project(v)
    print_output("Done updating tidal field")

solver_obj.load_state(start_time_point,outputdir='../../../outputs/0.validation/discrete-4cores')

cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')

#Effected area location
E_area_centre_point1 = [443385,3323260]
E_area_circle1 = 60

E_area_centre_point2 = [443786,3322300]
E_area_circle2 = 60

# Operation of tidal turbine farm about each turbine output through a callback
cb3 = rmse_r2.RMSECallback(solver_obj,'../../../outputs/0.validation/discrete-4cores', E_area_centre_point1, E_area_circle1)
solver_obj.add_callback(cb3,'timestep')

# Operation of tidal turbine farm about each turbine output through a callback
cb4 = rmse_r3.RMSECallback(solver_obj,'../../../outputs/0.validation/discrete-4cores', E_area_centre_point2, E_area_circle2)
solver_obj.add_callback(cb4,'timestep')

solver_obj.iterate(update_forcings=update_forcings)
power_output= sum(cb.average_power)
effect_two = (cb3.RMSEaverage + cb4.RMSEaverage)/2

###cable length###
turbine_locations = [float(x) for xy in farm_options.turbine_coordinates for x in xy]
landpointlocation = [444000,3323000]
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank == 0:
    cableclass = prepare_cable.Hybrid_Code.CableCostGA(turbine_locations=turbine_locations, substation_location=landpointlocation,capacity = 4)
    order_w = cableclass.compute_cable_cost_order()
else:
    order_w = []
order_w = comm.bcast(order_w, root=0)

landpointlocation_con = [Constant(x) for x in landpointlocation]
order_con = [Constant(i) for j in order_w for i in j]
cablecost = cablelength([x for xy in farm_options.turbine_coordinates for x in xy],landpointlocation_con,order_con)

shortestline = sqrt((444000-(443592-20))**2 + (3323000-(3322848-20))**2)


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank == 0:
    with open('result-l_y.txt','a+') as f:
        f.write(str(BEcable)+'\t'+str(BEsediment)+'\t')
        f.write(str(cb.average_power[-1])+'\t'+str(effect_two)+'\t'+str(cablecost)+'\t'+str(shortestline)+'\t'+str(iteration_numbers)+'\n')
        # f.write(str(iteration_numbers) +'\n')
else:
    pass

# print_output(str(BE)+'\t'+str(BE_sediment)+'\t')
# print_output(str(scaled_functional)+'\t'+str((cb.average_power[-1]-cb.average_profit[-1])/BE*10)+'\t'+str(cb.average_power[-1])+'\t'+str(effect_two)+'\t')


end_time = time.time()
print_output('time cost: {0:.2f}h'.format((end_time - start_time)/60/60))

if 1:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0 :
        yag = yagmail.SMTP(user = '623001493@qq.com',password = 'tiypryoerjyhbeib', host = 'smtp.qq.com')
        yag.send(to = ['623001493@qq.com'], subject = output_dir[17:], contents = ['Time cose: {0:.2f}h.'.format((end_time-start_time)/60/60)])
    else:
        pass

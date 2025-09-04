from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
import numpy
op2.init(log_level=INFO)
import sys
sys.path.append('../..')
import prepare_continuous.utm, prepare_continuous.myboundary, prepare_continuous.detectors
import os
import time
import yagmail
import rmse_sediment

start_time = time.time()

get_index = os.path.basename(sys.argv[0])
break_point = get_index.index('-')
BE = float(get_index[break_point+1:-3])
BE_sediment = float(get_index[:break_point])
print(str(BE)[:-2]+'_'+str(BE_sediment)[:-2])
# BE,BE_sediment = 45.0,0.0

###spring:676,middle:492,neap:340###
start_time_point = 486

output_dir = '../../../outputs/3.environment/zhoushan-continuous-forward/behind_notfrom0/'+str(BE)[:-2]+'_'+str(BE_sediment)[:-2]
print_output(output_dir[17:])



file_dir = '../../'
mesh2d = Mesh(file_dir+'mesh_continuous/mesh.msh')
#timestepping options
dt = 5*60 # reduce this if solver does not converge
t_export = 30*60 
# t_end = 1555200
t_end = 492*t_export+ 13*60*60


test_gradient, optimise = 0,1

turbine_area_PhyID = 2

P1 = FunctionSpace(mesh2d, "CG", 1)

# read bathymetry code
chk = DumbCheckpoint(file_dir+'prepare_continuous/bathymetry', mode=FILE_READ)
bathymetry2d = Function(P1)
chk.load(bathymetry2d, name='bathymetry')
chk.close()

#read viscosity / manning boundaries code
chk = DumbCheckpoint(file_dir+'prepare_continuous/viscosity', mode=FILE_READ)
h_viscosity = Function(P1, name='viscosity')
chk.load(h_viscosity)
chk.close()

# #manning = Function(P1,name='manning')
# chk = DumbCheckpoint(file_dir+'prepare_continuous/manning', mode=FILE_READ)
# manning = Function(bathymetry2d.function_space(), name='manning')
# chk.load(manning)
# chk.close()

#account for Coriolis code
def coriolis(mesh, lat,):
    R = 6371e3
    Omega = 7.292e-5
    lat_r = lat * pi / 180.
    f0 = 2 * Omega * sin(lat_r)
    beta = (1 / R) * 2 * Omega * cos(lat_r)
    x = SpatialCoordinate(mesh)
    x_0, y_0, utm_zone, zone_letter = prepare_continuous.utm.from_latlon(lat, 0)
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



# a density function (to be optimised below) that specifies the number of turbines per unit area
turbine_density = Function(get_functionspace(mesh2d, "CG", 1), name='turbine_density')
# associate subdomain_id 2 (as in dx(2)) with a tidal turbine farm
# (implemented via a drag term) with specified turbine density
# Turbine characteristic can be specified via:
# - farm_options.turbine_options.thrust_coefficient (default 0.8)
# - farm_options.turbine_options.diameter (default 16.0)
farm_options = TidalTurbineFarmOptions()
farm_options.turbine_density = turbine_density
# amount of power produced per turbine (kW) on average to "break even" (cost = revenue)
# this is used to scale the cost, which is assumed to be linear with the number of turbines,
# in such a way that the cost is expressed in kW which can be subtracted from the profit
# which is calculated as the power extracted by the turbines
farm_options.break_even_wattage = BE/10
options.tidal_turbine_farms[turbine_area_PhyID] = farm_options

# we first run the "forward" model with no turbines
# turbine_density.assign(0.0)
chk = DumbCheckpoint('../../../outputs/3.environment/zhoushan-continuous-op/behind_notfrom0/'+str(BE)[:-2]+'_'+str(BE_sediment)[:-2]+'/optimal_density', mode=FILE_READ)
# chk = DumbCheckpoint('../../../outputs/2.economy/continuous/intermediate/BE'+str(BE)[:-2]+'/optimal_density', mode=FILE_READ)
op_t_d = Function(bathymetry2d.function_space(), name='optimal_density')
chk.load(op_t_d)
chk.close()
turbine_density.project(op_t_d)

# create a density restricted to the farm
# the turbine_density, which is the control that will be varied in the optimisation,
# is defined everywhere, but it's influence is restricted to the farm area -
# the turbine drag term is integrated over the farm area only
# Because the turbine density is CG, the nodal values at the farm boundaries
# introduces a jagged edge around the farm where these values are tapered to zero.
# These nonzero values outside the farm itself do not contribute in any of the
# computations. For visualisation purposes we therefore project to a DG field
# restricted to the farm.
farm_density = Function(get_functionspace(mesh2d, "DG", 1), name='farm_density')
projector = SubdomainProjector(turbine_density, farm_density, turbine_area_PhyID)
projector.project()

cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')

# set boundary/initial conditions code
def update_forcings(t):
  with timed_stage('update forcings'):
    print_output("Updating tidal field at t={}".format(t))
    elev = prepare_continuous.myboundary.set_tidal_field(Function(bathymetry2d.function_space()), t, dt)
    tidal_elev.project(elev) 
    v = prepare_continuous.myboundary.set_velocity_field(Function(VectorFunctionSpace(mesh2d,"CG",1)),t,dt)
    tidal_v.project(v)
    print_output("Done updating tidal field")


tidal_elev = Function(bathymetry2d.function_space())
tidal_v = Function(VectorFunctionSpace(mesh2d,"CG",1))
solver_obj.bnd_functions['shallow_water'] = {
        200: {'elev': tidal_elev,'uv':tidal_v},  #set open boundaries to tidal_elev function
  }

solver_obj.load_state(start_time_point,outputdir='../../../outputs/0.validation/continuous-4cores')

# #place detectors code
# locations, names = prepare_continuous.detectors.get_detectors(mesh2d)
# cb_d = DetectorsCallback(solver_obj, locations, ['elev_2d', 'uv_2d','sediment_2d'], name='detectors',detector_names=names)
# solver_obj.add_callback(cb_d, 'timestep')

#Effected area location
E_area_centre_point = [[444500,3320800]] #behind
E_area_circle = [200]

# Operation of tidal turbine farm about each turbine output through a callback
cb2 = rmse_sediment.RMSECallback(solver_obj,'../../../outputs/3.environment/zhoushan-continuous/486-d5-4', E_area_centre_point, E_area_circle,492)
solver_obj.add_callback(cb2,eval_interval='export')

# #Effected area location
# E_area_centre_point = [[442700,3323760]] #inside
# E_area_circle = [100]

# # Operation of tidal turbine farm about each turbine output through a callback
# cb3 = rmse_sediment.RMSECallback(solver_obj,'../../../outputs/3.environment/zhoushan-continuous/486-d5-4', E_area_centre_point, E_area_circle,492)
# solver_obj.add_callback(cb3,eval_interval='export')

# #Effected area location
# E_area_centre_point1 = [[442973,3321135]] #next to
# E_area_circle1 = [200]

# # Operation of tidal turbine farm about each turbine output through a callback
# cb4 = rmse_sediment.RMSECallback(solver_obj,'../../../outputs/3.environment/zhoushan-continuous/486-d5-4', E_area_centre_point1, E_area_circle1,492)
# solver_obj.add_callback(cb4,eval_interval='export')

solver_obj.iterate(update_forcings=update_forcings)


# # compute maximum turbine density (based on a minimum of 2D between
# # turbines laterally, and 5D in the streamwise direction)
# D = farm_options.turbine_options.diameter
# max_density = 1./(2*D*2*D)
# print_output("Maximum turbine density = {}".format(max_density))

# we rescale the functional such that the gradients are ~ order magnitude 1.
# the scaling is based on the maximum cost term
# also we multiply by -1 so that if we minimize the functional, we maximize profit
# (maximize is also availble from pyadjoint but currently broken)
# scaling = -1/assemble(max(farm_options.break_even_wattage, 100) * max_density * dx(turbine_area_PhyID, domain=mesh2d))
# scaled_functional = scaling * cb.integrated_power
effect_two = [cb2.RMSEaverage]
if BE - 5 == 0:
    profit_max,sediment_max = 30979.69821,  54513.03909
elif BE - 25 == 0:
    profit_max,sediment_max = 29129.82268,	53619.17313
elif BE - 45 == 0:
    profit_max,sediment_max = 27355.68409,	51576.50981
scaled_functional = (1-BE_sediment/100)*cb.average_profit[-1]/profit_max*sediment_max- effect_two[0] * BE_sediment/100
print_output(effect_two)
print_output(cb.average_profit[-1])
print_output(scaled_functional)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank == 0:
    with open('result-l_y.txt','a+') as f:
        f.write(str(BE)+'\t'+str(BE_sediment)+'\t')
        f.write(str((cb.average_power[-1]-cb.average_profit[-1])/BE*10)+'\t'+str(cb.average_power[-1])+'\t'+str(effect_two[0])+'\n')
        # f.write(str(iteration_numbers) +'\n')
else:
    pass

print_output(str(BE)+'\t'+str(BE_sediment)+'\t')
print_output(str(scaled_functional)+'\t'+str((cb.average_power[-1]-cb.average_profit[-1])/BE*10)+'\t'+str(cb.average_power[-1])+'\t'+str(effect_two)+'\t')


end_time = time.time()
print_output('time cost: {0:.2f}h'.format((end_time - start_time)/60/60))

if 0:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0 :
        yag = yagmail.SMTP(user = '623001493@qq.com',password = 'tiypryoerjyhbeib', host = 'smtp.qq.com')
        yag.send(to = ['623001493@qq.com'], subject = output_dir[17:], contents = ['Time cose: {0:.2f}h.'.format((end_time-start_time)/60/60)])
    else:
        pass

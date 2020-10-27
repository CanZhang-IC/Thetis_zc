from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
import numpy
op2.init(log_level=INFO)
import sys
sys.path.append('..')
import prepare.utm, prepare.myboundary_30min
import time

t_start = time.time()

output_dir = '../../outputs/eachoutput_optimised_middle_30min-all'

mesh2d = Mesh('../mesh/mesh.msh')
#timestepping options
dt = 30*60 # reduce this if solver does not converge
t_export = 30*60 
#t_end = 1555200
#t_end = 1216800+ 13*60*60 # spring
t_end = 885600 + 13*60*60 # middle
#t_end = 612000 + 13*60*60 # neap
#t_end = 30*60


P1 = FunctionSpace(mesh2d, "CG", 1)

# read bathymetry code
chk = DumbCheckpoint('../prepare/bathymetry', mode=FILE_READ)
bathymetry2d = Function(P1)
chk.load(bathymetry2d, name='bathymetry')
chk.close()

#read viscosity / manning boundaries code
chk = DumbCheckpoint('../prepare/viscosity', mode=FILE_READ)
h_viscosity = Function(P1, name='viscosity')
chk.load(h_viscosity)
chk.close()

#manning = Function(P1,name='manning')
chk = DumbCheckpoint('../prepare/manning', mode=FILE_READ)
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
        200: {'elev': tidal_elev},  #set open boundaries to tidal_elev function
  }

def update_forcings(t):
  with timed_stage('update forcings'):
    print_output("Updating tidal field at t={}".format(t))
    elev = prepare.myboundary_30min.set_tidal_field(Function(bathymetry2d.function_space()), t, dt)
    tidal_elev.project(elev) 
    v = prepare.myboundary_30min.set_velocity_field(Function(VectorFunctionSpace(mesh2d,"CG",1)),t,dt)
    tidal_v.project(v)
    print_output("Done updating tidal field")


# Initialise Discrete turbine farm characteristics
farm_options = DiscreteTidalTurbineFarmOptions()
farm_options.turbine_type = 'constant'
farm_options.turbine_options.thrust_coefficient = 0.6
farm_options.turbine_options.diameter = 20
farm_options.upwind_correction = False

site_x1, site_y1, site_x2, site_y2 = 443342 ,3322632, 443591, 3322845


result_data = [
    #coordinates
    443352.0014672894, 3322661.2249424085, 443352.00360414607, 3322799.4023007597, 443370.24033160694, 3322835.0, 443412.747696928, 3322834.998605353, 443373.21377413854, 3322765.489792042, 443391.9741800234, 3322800.817274995, 443433.9590554587, 3322801.0884508095, 443455.17338239006, 3322834.9980192194, 443413.1878277575, 3322766.9040195737, 443456.8741346038, 3322768.30120755, 443478.0871860666, 3322802.2139794477, 443500.9999548514, 3322834.999785881, 443501.1198651514, 3322769.510135885, 443530.80403638876, 3322796.3215907435, 443540.9999975333, 3322834.998255678, 443580.99752880324, 3322834.9969152194, 
    #axis
    99.31049255675143, 97.79287314524697, 100.21135081450862, 101.81294766091462, 98.1446667237933, 99.56009623162055, 101.43650213959047, 103.91391264074014, 100.02038449487331, 101.9470124097927, 103.31352373072824, 106.09853933715259, 103.42342515732167, 105.4148131263235, 107.67334339092514, 109.17335840499666
    ]


farm_options.turbine_coordinates = [[Constant(result_data[2*i]), Constant(result_data[2*i+1])] for i in range(int(len(result_data)/3))]
farm_options.considering_yaw = True
farm_options.turbine_axis = [Constant(i) for i in result_data[int(len(result_data)/3*2):]]
#add turbines to SW_equations
options.discrete_tidal_turbine_farms[2] = farm_options

###spring:676,middle:492,neap:340###
solver_obj.load_state(492, outputdir='../../outputs/redata_30min_normaldepth')
#solver_obj.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=Constant(0.0))

# Operation of tidal turbine farm through a callback
cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')

cb2 = turbines.EachTurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb2, 'timestep')
# start computer forward model

solver_obj.iterate(update_forcings=update_forcings)

power_output= sum(cb.integrated_power)
print(power_output)

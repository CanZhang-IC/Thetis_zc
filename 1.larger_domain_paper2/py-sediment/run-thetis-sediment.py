from thetis import *
import detectors
import tidal_forcing
import utm
import yagmail
import time

time_start= time.time()

output_dir = '../outputs/5min-4cores-sediment-exner'
mesh2d = Mesh('../mesh/mesh.msh')
#timestepping options
dt = 5*60 # reduce this if solver does not converge
t_export = 5*60 
t_end = 1555200+86400 # e.g. 16days+ 2day spin up = 1382400 s + 172800s = 1555200 s
#640800: 16/08/2013 09:59

P1 = FunctionSpace(mesh2d, "CG", 1)

# read bathymetry code
chk = DumbCheckpoint('bathymetry', mode=FILE_READ)
bathymetry2d = Function(P1)
chk.load(bathymetry2d, name='bathymetry')
chk.close()

#read viscosity / manning boundaries code
chk = DumbCheckpoint('viscosity', mode=FILE_READ)
h_viscosity = Function(P1, name='viscosity')
chk.load(h_viscosity)
chk.close()

#manning = Constant(0.02)
# chk = DumbCheckpoint('manning', mode=FILE_READ)
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
    x_0, y_0, utm_zone, zone_letter = utm.from_latlon(lat, 0)
    coriolis_2d = Function(FunctionSpace(mesh, 'CG', 1), name="coriolis_2d")
    coriolis_2d.interpolate(f0 + beta * (x[1] - y_0))
    return coriolis_2d
coriolis_2d =coriolis(mesh2d, 29)


# --- create solver ---
solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry2d)
options = solver_obj.options

sediment_size,morfac,diffusivity,morphological_viscosity = 1e-3,1,1,1e-9
options.sediment_model_options.solve_suspended_sediment = True
options.sediment_model_options.use_bedload = True
options.sediment_model_options.solve_exner = True
# options.sediment_model_options.use_angle_correction = True
# options.sediment_model_options.use_slope_mag_correction = True
# options.sediment_model_options.use_secondary_current = True
# options.sediment_model_options.use_advective_velocity_correction = False

options.sediment_model_options.use_sediment_conservative_form = False
options.sediment_model_options.average_sediment_size = Constant(sediment_size)
options.sediment_model_options.bed_reference_height = Constant(0.025)
options.sediment_model_options.morphological_acceleration_factor = Constant(morfac)
options.sediment_model_options.morphological_viscosity = Constant(morphological_viscosity)
options.sediment_model_options.porosity = Constant(0.4)
options.horizontal_viscosity = h_viscosity
# define parameters
options.nikuradse_bed_roughness = Constant(3*sediment_size)
options.horizontal_diffusivity = Constant(diffusivity)
options.norm_smoother = Constant(0.1)

options.cfl_2d = 1.0
options.use_nonlinear_equations = True
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.coriolis_frequency = coriolis_2d
options.output_directory =  output_dir
options.check_volume_conservation_2d = True
if options.sediment_model_options.solve_suspended_sediment:
    options.fields_to_export = ['sediment_2d', 'uv_2d', 'elev_2d', 'bathymetry_2d']  # note exporting bathymetry must be done through export func
    options.fields_to_export_hdf5 = ['sediment_2d', 'uv_2d', 'elev_2d', 'bathymetry_2d']
    options.sediment_model_options.check_sediment_conservation = True
else:
    options.fields_to_export = ['uv_2d', 'elev_2d', 'bathymetry_2d']  # note exporting bathymetry must be done through export func
    options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d', 'bathymetry_2d']
options.element_family = "dg-dg"
options.timestepper_type = 'CrankNicolson'
options.timestepper_options.implicitness_theta = 1 #Implicitness parameter, default 0.5
options.timestepper_options.use_semi_implicit_linearization = True#If True use a linearized semi-implicit scheme
options.use_wetting_and_drying = True #Wetting and drying is included through the modified bathymetry formulation of Karna et al. (2011). 
options.wetting_and_drying_alpha = Constant(0.5) #need to check if this is a good value
# options.manning_drag_coefficient = manning 
#options.quadratic_drag_coefficient = Constant(0.0015)

#Epshteyn and Riviere (2007). Estimation of penalty parameters for symmetric interior penalty Galerkin methods. Journal of Computational and Applied Mathematics, 206(2):843-872. http://dx.doi.org/10.1016/j.cam.2006.08.029
# options.use_grad_div_viscosity_term = True
# options.use_grad_depth_viscosity_term = False

options.timestep = dt

options.timestepper_options.solver_parameters = {'snes_monitor': None,
                                                 'snes_rtol': 1e-9,
                                                 'ksp_type': 'preonly',
                                                 'pc_type': 'lu',
                                                 'pc_factor_mat_solver_type': 'mumps',
                                                 'mat_type': 'aij'
                                                 }


# set boundary/initial conditions code
tidal_elev = Function(bathymetry2d.function_space())
solver_obj.bnd_functions['shallow_water'] = {
        200: {'elev': tidal_elev},  #set open boundaries to tidal_elev function
  }

# initialise velocity and elevation
DG_2d = FunctionSpace(mesh2d, "DG", 1)
vector_dg = VectorFunctionSpace(mesh2d, "DG", 1)

chk = DumbCheckpoint("../outputs/5min-4cores-sediment-hydro/hdf5/Elevation2d_00020", mode=FILE_READ)
elev = Function(DG_2d, name="elev_2d")
chk.load(elev)
chk.close()

chk = DumbCheckpoint("../outputs/5min-4cores-sediment-hydro/hdf5/Velocity2d_00020", mode=FILE_READ)
uv = Function(vector_dg, name="uv_2d")
chk.load(uv)
chk.close()

if options.sediment_model_options.solve_suspended_sediment:
    # setting an equilibrium boundary conditions results in the sediment value at the boundary
    # being chosen so that erosion and deposition are equal here (ie. in equilibrium) and the bed is immobile at this boundary
    solver_obj.bnd_functions['sediment'] = {
        200: {'elev': tidal_elev, 'equilibrium': None}}
    # set initial conditions
    solver_obj.assign_initial_conditions(uv=uv, elev=elev)

else:
    # set initial conditions
    solver_obj.assign_initial_conditions(uv=uv, elev=elev)
# solver_obj.assign_initial_conditions(uv=uv,elev=elev)
# restarting
# solver_obj.load_state(288,outputdir='../outputs/5min-1core-sediment-hydro')

#place detectors code
locations, names = detectors.get_detectors(mesh2d)
cb = DetectorsCallback(solver_obj, locations, ['elev_2d', 'uv_2d','sediment_2d'], name='detectors',detector_names=names)
solver_obj.add_callback(cb, 'timestep')

def update_forcings(t):
    with timed_stage('update forcings'):
        # print_output("Updating tidal field at t={}".format(t))
        tidal_forcing.set_tidal_field(tidal_elev, t)
        # print_output("Done updating tidal field")


solver_obj.iterate(update_forcings=update_forcings)

time_end= time.time()


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank == 0:
    yag = yagmail.SMTP(user = '623001493@qq.com',password = 'ouehigyjxpidbbcj', host = 'smtp.qq.com')
    yag.send(to = ['623001493@qq.com'], subject = 'Python done', contents = ['College computer Python Finished'])
    print('Time cost: {0:.2f}h'.format((time_end-time_start)/60/60))
else:
    pass
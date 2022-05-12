from thetis import *
from firedrake_adjoint import *
# op2.init(log_level=INFO)
import sys
sys.path.append('..')
import prepare.utm, prepare.myboundary, prepare.detectors
import time
import yagmail

start_time = time.time()

ouput_dir = './outputs/validation_5min-e&v-sediment'

mesh2d = Mesh('./mesh2/mesh.msh')
#timestepping options
dt = 5*60 # reduce this if solver does not converge
t_export = 5*60 
t_end = 1350000#1555200
# t_end = 1216800+ 13*60*60 # spring
# t_end = 885600 + 13*60*60 # middle
#t_end = 612000 + 13*60*60 # neap
#t_end = 30*60


P1 = FunctionSpace(mesh2d, "CG", 1)

# read bathymetry code
chk = DumbCheckpoint('./prepare/bathymetry', mode=FILE_READ)
bathymetry2d = Function(P1)
chk.load(bathymetry2d, name='bathymetry')
chk.close()

#read viscosity / manning boundaries code
chk = DumbCheckpoint('./prepare/viscosity', mode=FILE_READ)
h_viscosity = Function(P1, name='viscosity')
chk.load(h_viscosity)
chk.close()

#manning = Function(P1,name='manning')
# chk = DumbCheckpoint('./prepare/manning', mode=FILE_READ)
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
    x_0, y_0, utm_zone, zone_letter = prepare.utm.from_latlon(lat, 0)
    coriolis_2d = Function(FunctionSpace(mesh, 'CG', 1), name="coriolis_2d")
    coriolis_2d.interpolate(f0 + beta * (x[1] - y_0))
    return coriolis_2d
coriolis_2d =coriolis(mesh2d, 30)



# --- create solver ---
solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry2d)
options = solver_obj.options

sediment_size,morfac,diffusivity,morphological_viscosity = 5*1e-5,1,0.15,1e-6
options.sediment_model_options.solve_suspended_sediment = True
# options.sediment_model_options.use_bedload = True
options.sediment_model_options.solve_exner = True
# options.sediment_model_options.use_angle_correction = True
# options.sediment_model_options.use_slope_mag_correction = True
# options.sediment_model_options.use_secondary_current = True
# options.sediment_model_options.use_advective_velocity_correction = False

options.sediment_model_options.use_sediment_conservative_form = True
options.sediment_model_options.average_sediment_size = Constant(sediment_size)
options.sediment_model_options.bed_reference_height = Constant(3*sediment_size)
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
options.output_directory = ouput_dir
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
#options.quadratic_drag_coefficient = Constant(0.0025)

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
        1000: {'elev': tidal_elev,'uv':tidal_v},  #set open boundaries to tidal_elev function
  }

# initialise velocity and elevation
DG_2d = FunctionSpace(mesh2d, "DG", 1)
vector_dg = VectorFunctionSpace(mesh2d, "DG", 1)

chk = DumbCheckpoint("./outputs/validation_5min-e&v-sediment_hydro/hdf5/Elevation2d_00024", mode=FILE_READ)
elev = Function(DG_2d, name="elev_2d")
chk.load(elev)
chk.close()

chk = DumbCheckpoint("./outputs/validation_5min-e&v-sediment_hydro/hdf5/Velocity2d_00024", mode=FILE_READ)
uv = Function(vector_dg, name="uv_2d")
chk.load(uv)
chk.close()

if options.sediment_model_options.solve_suspended_sediment:
    # setting an equilibrium boundary conditions results in the sediment value at the boundary
    # being chosen so that erosion and deposition are equal here (ie. in equilibrium) and the bed is immobile at this boundary
    solver_obj.bnd_functions['sediment'] = {
        1000: {'elev': tidal_elev,'uv':tidal_v, 'equilibrium': None}}
    # set initial conditions
    solver_obj.assign_initial_conditions(uv=uv, elev=elev)

else:
    # set initial conditions
    solver_obj.assign_initial_conditions(uv=uv, elev=elev)
    pass

###spring:676,middle:492,neap:340###
# solver_obj.load_state(606, outputdir='./outputs/validation_5min-e&v-1')
# solver_obj.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=Constant(0.0))

#place detectors code
with stop_annotating():
  locations, names = prepare.detectors.get_detectors(mesh2d)
cb = DetectorsCallback(solver_obj, locations, ['elev_2d', 'uv_2d'], name='detectors',detector_names=names)
solver_obj.add_callback(cb, 'timestep')

def update_forcings(t):
  with timed_stage('update forcings'):
    # print_output("Updating tidal field at t={}".format(t))
    elev = prepare.myboundary.set_tidal_field(Function(bathymetry2d.function_space()), t, dt)
    tidal_elev.project(elev) 
    v = prepare.myboundary.set_velocity_field(Function(VectorFunctionSpace(mesh2d,"CG",1)),t,dt)
    tidal_v.project(v)
    # print_output("Done updating tidal field")

# start computer forward model
solver_obj.iterate(update_forcings=update_forcings)

end_time = time.time()

print('The time cost is {0:.2f}min'.format((end_time-start_time)/60))

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank == 0:
    yag = yagmail.SMTP(user = '623001493@qq.com',password = 'ouehigyjxpidbbcj', host = 'smtp.qq.com')
    yag.send(to = ['623001493@qq.com'], subject = 'Python done', contents = ['College computer Python Finished'])
else:
    pass
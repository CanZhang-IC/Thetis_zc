from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
import numpy
op2.init(log_level=INFO)
import sys
sys.path.append('../..')
import prepare_continuous.utm, prepare_continuous.tidal_forcing, prepare_continuous.detectors
import os
import time
import yagmail

t_start = time.time()

file_dir = '../../'

sediment_size,morfac,diffusivity,morphological_viscosity = 40*1e-6,1,0.01,1e-6
output_dir = '../../../outputs/0.validation/continuous-sediment-exner-tideforcing'
print_output(output_dir[17:])

mesh2d = Mesh(file_dir+'mesh_continuous/mesh.msh')

#timestepping options
dt = 5*60 # reduce this if solver does not converge
t_export = 30*60 
t_end = 1555200
# t_end = 1216800+ 13*60*60 # spring
# t_end = 885600 + 13*60*60 # middle
# t_end = 612000 + 13*60*60 # neap

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

options.sediment_model_options.solve_suspended_sediment = True
options.sediment_model_options.use_bedload = False
options.sediment_model_options.solve_exner = False
# options.sediment_model_options.use_angle_correction = True
# options.sediment_model_options.use_slope_mag_correction = True
# options.sediment_model_options.use_secondary_current = True
# options.sediment_model_options.use_advective_velocity_correction = False

options.sediment_model_options.use_sediment_conservative_form = True
options.sediment_model_options.average_sediment_size = sediment_size
options.sediment_model_options.bed_reference_height = 3*sediment_size
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
    options.fields_to_export = ['sediment_2d', 'uv_2d', 'elev_2d', 'bathymetry_2d']  # note exporting bathymetry must be done through export func
    options.fields_to_export_hdf5 = ['sediment_2d', 'uv_2d', 'elev_2d']
    options.sediment_model_options.check_sediment_conservation = True
else:
    options.fields_to_export = ['uv_2d', 'elev_2d', 'bathymetry_2d']  # note exporting bathymetry must be done through export func
    options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']
options.element_family = "dg-dg"
options.timestepper_type = 'CrankNicolson'
options.timestepper_options.implicitness_theta = 1 #Implicitness parameter, default 0.5
options.timestepper_options.use_semi_implicit_linearization = True#If True use a linearized semi-implicit scheme
options.use_wetting_and_drying = True #Wetting and drying is included through the modified bathymetry formulation of Karna et al. (2011). 
options.wetting_and_drying_alpha = Constant(0.5) #need to check if this is a good value
# options.manning_drag_coefficient = manning 
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
    prepare_continuous.tidal_forcing.set_tidal_field(tidal_elev, t)
    print_output("Done updating tidal field")

# # initialise velocity and elevation
# DG_2d = FunctionSpace(mesh2d, "DG", 1)
# vector_dg = VectorFunctionSpace(mesh2d, "DG", 1)

# chk = DumbCheckpoint("../../../outputs/1.larger_domain_paper2/sediment-hydro/hdf5/Elevation2d_00048", mode=FILE_READ)
# elev = Function(DG_2d, name="elev_2d")
# chk.load(elev)
# chk.close()

# chk = DumbCheckpoint("../../../outputs/1.larger_domain_paper2/sediment-hydro/hdf5/Velocity2d_00048", mode=FILE_READ)
# uv = Function(vector_dg, name="uv_2d")
# chk.load(uv)
# chk.close()

# if options.sediment_model_options.solve_suspended_sediment:
#     # setting an equilibrium boundary conditions results in the sediment value at the boundary
#     # being chosen so that erosion and deposition are equal here (ie. in equilibrium) and the bed is immobile at this boundary
#     solver_obj.bnd_functions['sediment'] = {
#         200: {'elev': tidal_elev, 'equilibrium': None}}
#     # set initial conditions
#     solver_obj.assign_initial_conditions(uv=uv, elev=elev)

# else:
#     # set initial conditions
#     solver_obj.assign_initial_conditions(uv=uv, elev=elev)
###spring:676,middle:492,neap:340###
solver_obj.assign_initial_conditions(uv=as_vector((1.0, 0.0)), elev=Constant(0.5))
# solver_obj.load_state(574,outputdir='../../../outputs/0.validation/continuous-sediment-exner')

#place detectors code
with stop_annotating():
  locations, names = prepare_continuous.detectors.get_detectors(mesh2d)
cb = DetectorsCallback(solver_obj, locations, ['elev_2d', 'uv_2d','sediment_2d'], name='detectors',detector_names=names)
solver_obj.add_callback(cb, 'timestep')
# start computer forward model
solver_obj.iterate(update_forcings=update_forcings)

t_end = time.time()

print('The time cost is {0:.2f} min'.format((t_end-t_start)/60))

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank == 0:
    yag = yagmail.SMTP(user = '623001493@qq.com',password = 'ouehigyjxpidbbcj', host = 'smtp.qq.com')
    yag.send(to = ['623001493@qq.com'], subject = 'My computer', contents = [output_dir+' ###### '+ 'Time cose: {0:.2f}h.'.format((t_end-t_start)/60/60)])
else:
    pass
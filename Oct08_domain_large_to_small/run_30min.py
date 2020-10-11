from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
op2.init(log_level=INFO)
import detectors
import tidal_forcing
import utm
import scipy.interpolate
import myboundary_30min
import time


time_start=time.time()

ouput_dir = 'outputs/redata-dgdg-Latitude-30'


mesh2d = Mesh('./mesh/mesh.msh')
#timestepping options
dt = 30*60 # reduce this if solver does not converge
t_export = 30*60 
t_end = 1555200#288*60*60


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

#manning = Function(P1,name='manning')
chk = DumbCheckpoint('manning', mode=FILE_READ)
manning = Function(bathymetry2d.function_space(), name='manning')
chk.load(manning)
chk.close()
# mesh_manning = Mesh('../zs_manning/mesh/mesh.msh')
# chk = DumbCheckpoint('optimal_manning_neap', mode=FILE_READ)
# manning_read = Function(FunctionSpace(mesh_manning,"CG",1), name='optimal_manning_neap')
# chk.load(manning_read)
# chk.close()

# xvector1 = mesh_manning.coordinates.dat.data
# mvector1 = manning_read.dat.data
# assert xvector1.shape[0] == mvector1.shape[0]
# xx=[]
# yy=[]
# for xy in xvector1:
#   xx.append(xy[0])
#   yy.append(xy[1])

# manning = Function(bathymetry2d.function_space(), name='manning')
# xvector2 = mesh2d.coordinates.dat.data
# interpolator = scipy.interpolate.griddata(xvector1,mvector1,xvector2)
# manning.project(manning_read)
# mvector2 = manning.dat.data
# assert xvector2.shape[0] == mvector2.shape[0]
# for i,xy in enumerate(xvector2):
#   mvector2[i] = interpolator((xy[0],xy[1]))

# def smoothen_manning(manning2d):  # smoothing manning
#   v = TestFunction(manning2d.function_space())

#   massb = assemble(v * manning2d *dx)#装配起来
#   massl = assemble(v*dx)
#   with massl.dat.vec as ml, massb.dat.vec as mb, manning2d.dat.vec as sb:
#       ml.reciprocal()
#       sb.pointwiseMult(ml, mb)

# smoothen_manning(manning)
# smoothen_manning(manning)
# File('read_manning.pvd').write(manning)


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
coriolis_2d =coriolis(mesh2d, 30)

# --- create solver ---
solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry2d)
options = solver_obj.options
options.cfl_2d = 1.0
options.use_nonlinear_equations = True
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.coriolis_frequency = coriolis_2d
options.output_directory = ouput_dir
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
        200: {'elev': tidal_elev, 'uv':tidal_v},  #set open boundaries to tidal_elev function
  }

def update_forcings(t):
  with timed_stage('update forcings'):
    print_output("Updating tidal field at t={}".format(t))
    elev = myboundary_30min.set_tidal_field(Function(bathymetry2d.function_space()), t,dt)
    tidal_elev.project(elev) 
    v = myboundary_30min.set_velocity_field(Function(VectorFunctionSpace(mesh2d,"CG",1)),t,dt)
    tidal_v.project(v)
    print_output("Done updating tidal field")

# run as normal (this run will be annotated by firedrake_adjoint)
#solver_obj.load_state(492,outputdir='./outputs/test')
solver_obj.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=Constant(0.0))

#place detectors code
with stop_annotating():
  locations, names = detectors.get_detectors(mesh2d)
cb = DetectorsCallback(solver_obj, locations, ['elev_2d', 'uv_2d'], name='detectors',detector_names=names)
solver_obj.add_callback(cb, 'timestep')

solver_obj.iterate(update_forcings=update_forcings)

time_end=time.time()
print('time cost',time_end-time_start,'s')
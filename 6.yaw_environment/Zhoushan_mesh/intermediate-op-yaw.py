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

t_start = time.time()

file_dir = '../../'

# get_index = os.path.basename(sys.argv[0])
# namelength = len('paper3_run')
# P_factor = float(get_index[namelength:-3])

P_factor = 0.8
print_output(str(P_factor))
output_dir = '../../../outputs/6.yaw_environment/Paper3/Zhoushan_mesh/optimisation/backhome/intermediate-yaw_op-P_factor_'+str(P_factor)+'-5min_e&v'

mesh2d = Mesh(file_dir+'mesh/mesh.msh')

#timestepping options
dt = 5*60 # reduce this if solver does not converge
t_export = 30*60 
# t_end = 1555200
# t_end = 1216800+ 13*60*60 # spring
t_end = 885600 + 13*60*60 # middle
# t_end = 612000 + 13*60*60 # neap


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

turbine_location = []
for x in range(xmin+20,xmax-20,60):
    for y in range(ymin+20,ymax-20,120):
        turbine_location.append([x,y])
for x in range(xmin+20+30,xmax-20,60):
    for y in range(ymin+20+60,ymax-20,120):
        turbine_location.append([x,y])
farm_options.turbine_coordinates =[[Constant(xy[0]),Constant(xy[1])] for xy in turbine_location]

# farm_options.considering_yaw = True
# ### flood_dir,ebb_dir = [151.16498222165774, 158.5544801109472, 139.21937282396664, 138.08305879290157, 124.96583232750241, 115.81910043683706, 122.45039350653985, 90.7905468922789, 151.82514064875716, 137.42231283489394, 139.31962364099329, 107.92321470506059],[271.05995847732805, 358.4803421898285, 289.66268839560024, 290.57781969194, 301.2642266120518, 305.56975116938287, 334.1911304453266, 324.0491913809176, 450.0, 305.72560256885924, 315.6020927982799, 344.27854640341263]

# ### farm_options.turbine_axis = [Constant(i) for i in flood_dir] + [Constant(i) for i in ebb_dir]

# farm_options.turbine_axis = [Constant(180) for i in range(len(farm_options.turbine_coordinates))] + [Constant(360) for i in range(len(farm_options.turbine_coordinates))]

result_output_dir = '../../../outputs/6.yaw_environment/Paper3/Zhoushan_mesh/optimisation/backhome/intermediate-yaw_op-P_factor_'+str(0.7)+'-5min_e&v'
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

farm_options.considering_yaw = True
flood_dir,ebb_dir = all_controls[:12], all_controls[12:]

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

###spring:676,middle:492,neap:340###
# solver_obj.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=Constant(0.0))
solver_obj.load_state(492, outputdir='../../../outputs/6.yaw_environment/Paper3/Zhoushan_mesh/restart_5min-e&v')

# Operation of tidal turbine farm through a callback
cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')

#Effected area location
E_area_centre_point = [(xmin+xmax)/2,ymax+800]
E_area_circle = 60

# Operation of tidal turbine farm about each turbine output through a callback
cb2 = rmse_r2.RMSECallback(solver_obj,'../../../outputs/6.yaw_environment/Paper3/Zhoushan_mesh/optimisation/forward/intermediate-forward-5min_e&v', E_area_centre_point, E_area_circle)
solver_obj.add_callback(cb2,'timestep')

# start computer forward model
solver_obj.iterate(update_forcings=update_forcings)

# ###set up interest functional and control###
power_output= sum(cb.average_power)
maxoutput, maxeffect = 2306.956713596631,	3389.005152125812
interest_functional = (P_factor*(power_output/maxoutput)-(1-P_factor)*(cb2.RMSEaverage/maxeffect))*maxoutput
print(interest_functional,power_output,cb2.RMSEaverage)
# interest_functional = sum(cb.average_power)

# specifies the control we want to vary in the optimisation
c =[Control(x) for x in farm_options.turbine_axis]#[Control(x) for xy in farm_options.turbine_coordinates for x in xy]

# a number of callbacks to provide output during the optimisation iterations:
# - ControlsExportOptimisationCallback export the turbine_friction values (the control)
#            to outputs/control_turbine_friction.pvd. This can also be used to checkpoint
#            the optimisation by using the export_type='hdf5' option.
# - DerivativesExportOptimisationCallback export the derivative of the functional wrt
#            the control as computed by the adjoint to outputs/derivative_turbine_friction.pvd
# - UserExportOptimisationCallback can be used to output any further functions used in the
#            forward model. Note that only function states that contribute to the functional are
#            guaranteed to be updated when the model is replayed for different control values.
# - FunctionalOptimisationCallback simply writes out the (scaled) functional values
# - the TurbineOptimsationCallback outputs the average power, cost and profit (for each
#            farm if multiple are defined)
turbine_density = Function(solver_obj.function_spaces.P1_2d, name='turbine_density')
turbine_density.interpolate(solver_obj.tidal_farms[0].turbine_density)
callback_list = optimisation.OptimisationCallbackList([
    optimisation.ConstantControlOptimisationCallback(solver_obj, array_dim=len(c)),
    optimisation.DerivativeConstantControlOptimisationCallback(solver_obj, array_dim=len(c)),
    optimisation.UserExportOptimisationCallback(solver_obj, [turbine_density, solver_obj.fields.uv_2d]),
    optimisation.FunctionalOptimisationCallback(solver_obj),
])
# callbacks to indicate start of forward and adjoint runs in log
def eval_cb_pre(controls):
    print_output("FORWARD RUN:")
    print_output("Controls: {}".format([float(c) for c in controls]))

def derivative_cb_pre(controls):
    print_output("ADJOINT RUN:")
    print_output("Controls: {}".format([float(c) for c in controls]))

# this reduces the functional J(u, m) to a function purely of the control m:
# rf(m) = J(u(m), m) where the velocities u(m) of the entire simulation
# are computed by replaying the forward model for any provided turbine coordinates m
rf = ReducedFunctional(-interest_functional, c, derivative_cb_post=callback_list,
        eval_cb_pre=eval_cb_pre, derivative_cb_pre=derivative_cb_pre)

if 0:
    # whenever the forward model is changed - for example different terms in the equation,
    # different types of boundary conditions, etc. - it is a good idea to test whether the
    # gradient computed by the adjoint is still correct, as some steps in the model may
    # not have been annotated correctly. This can be done via the Taylor test.
    # Using the standard Taylor series, we should have (for a sufficiently smooth problem):
    #   rf(td0+h*dtd) - rf(td0) - < drf/dtd(rf0), h dtd> = O(h^2)

    # we choose a random point in the control space, i.e. a randomized turbine density with
    # values between 0 and 1 and choose a random direction dtd to vary it in

    # this tests whether the above Taylor series residual indeed converges to zero at 2nd order in h as h->0
   
    m0 =  [Constant(x) for x in farm_options.turbine_axis] # [Constant(x) for xy in farm_options.turbine_coordinates for x in xy] 
    h0 =  [Constant(1) for x in farm_options.turbine_axis] # [Constant(1) for xy in farm_options.turbine_coordinates for x in xy] 

    minconv = taylor_test(rf, m0, h0)
    print_output("Order of convergence with taylor test (should be 2) = {}".format(minconv))

    assert minconv > 1.95

if 1:
    # Optimise the control for minimal functional (i.e. maximum profit)
    # with a gradient based optimisation algorithm using the reduced functional
    # to replay the model, and computing its derivative via the adjoint
    # By default scipy's implementation of L-BFGS-B is used, see
    #   https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html
    # options, such as maxiter and pgtol can be passed on.
    optimise_layout_only = False
    d = farm_options.turbine_options.diameter

    lb = [90]*int(len(farm_options.turbine_axis)/2) + [270]*int(len(farm_options.turbine_axis)/2) 
    ub = [270]*int(len(farm_options.turbine_axis)/2) + [450]*int(len(farm_options.turbine_axis)/2) 

    
    td_opt = minimize(rf, method='SLSQP', bounds=[lb,ub],options={'maxiter': 100, 'ptol': 1e-3})


t_end = time.time()
print('time cost: {0:.2f}h'.format((t_end - t_start)/60/60))

if 1:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0 :
        yag = yagmail.SMTP(user = '623001493@qq.com',password = 'ouehigyjxpidbbcj', host = 'smtp.qq.com')
        yag.send(to = ['canzhang2019@gmail.com'], subject = 'Python done', contents = ['Intermediate yaw Optimisation cost: '+str((t_end - t_start)/60/60/24) + 'days.'+'P_factoe:'+str(P_factor)])
    else:
        pass
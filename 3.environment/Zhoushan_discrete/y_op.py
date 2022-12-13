from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
op2.init(log_level=INFO)
import sys
sys.path.append('../..')
import prepare.utm, prepare.myboundary
import prepare_cable.Hybrid_Code
import rmse_r2
import rmse_r3
import numpy
import time
import yagmail
import h5py

start_time = time.time()

get_index = os.path.basename(sys.argv[0])
BE = float(get_index[:-3])
# BE = 0.0


t1_end, t2_end = 503, 515
P_factor = 1.0
file_dir = '../../'
output_dir = '../../../outputs/3.environment/discrete/flood_ebb/y-BE-'+str(BE)[:-2]
print_output(output_dir[17:])
mesh2d = Mesh(file_dir+'mesh/mesh.msh')

#timestepping options
dt = 5*60 # reduce this if solver does not converge
t_export = 30*60 
# t_end = 1555200
# t_end = 1216800+ 13*60*60 # spring
t_end1 = t1_end*30*60 + 2*60*60# middle
t_end2 = t2_end*30*60 + 2*60*60
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


# Initialise Discrete turbine farm characteristics
farm_options = DiscreteTidalTurbineFarmOptions()
farm_options.turbine_type = 'constant'
farm_options.turbine_options.thrust_coefficient = 0.6
farm_options.turbine_options.diameter = 20
# farm_options.turbine_options.C_support = 1
# farm_options.turbine_options.A_support = H/2
farm_options.upwind_correction = True

xmin,ymin,xmax,ymax = 443340, 3322634, 443592, 3322848 


# turbine_location = []
# x_space = 60
# for x in range(xmin+20,xmax-20,x_space):
#     for y in range(ymin+20,ymax-20,x_space*2):
#         turbine_location.append([x,y])
# for x in range(xmin+20+int(x_space/2),xmax-20,x_space):
#     for y in range(ymin+20+x_space,ymax-20,x_space*2):
#         turbine_location.append([x,y])
# farm_options.turbine_coordinates =[[Constant(xy[0]),Constant(xy[1])] for xy in turbine_location]

# farm_options.considering_yaw = True
# farm_options.turbine_axis = [Constant(90) for i in range(len(farm_options.turbine_coordinates))] + [Constant(270) for i in range(len(farm_options.turbine_coordinates))]

result_output_dir = '../../../outputs/3.environment/discrete/flood_ebb/y-BE-'+str(BE-2)[:-2]
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
options2.discrete_tidal_turbine_farms[2] = options.discrete_tidal_turbine_farms[2]


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
solver_obj.load_state(t1_end, outputdir='../../../outputs/0.validation/discrete-4cores')
# Operation of tidal turbine farm through a callback
cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')

#Effected area location
E_area_centre_point1 = [442888,3323940]
E_area_circle1 = 60

E_area_centre_point2 = [443786,3322300]
E_area_circle2 = 60

# Operation of tidal turbine farm about each turbine output through a callback
cb3 = rmse_r2.RMSECallback(solver_obj,'../../../outputs/0.validation/discrete-4cores', E_area_centre_point1, E_area_circle1)
solver_obj.add_callback(cb3,'timestep')

# Operation of tidal turbine farm about each turbine output through a callback
cb4 = rmse_r3.RMSECallback(solver_obj,'../../../outputs/0.validation/discrete-4cores', E_area_centre_point2, E_area_circle2)
solver_obj.add_callback(cb4,'timestep')

# start computer forward model
solver_obj.iterate(update_forcings=update_forcings)
###set up interest functional and control###
power_output1= sum(cb.average_power)
effect_two1 = (cb3.RMSEaverage + cb4.RMSEaverage)/2

###spring:676,middle:492,neap:340###
# solver_obj.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=Constant(0.0))
solver_obj2.load_state(t2_end, outputdir='../../../outputs/0.validation/discrete-4cores')
# Operation of tidal turbine farm through a callback
cb2 = turbines.TurbineFunctionalCallback(solver_obj2)
solver_obj2.add_callback(cb2, 'timestep')

# Operation of tidal turbine farm about each turbine output through a callback
cb5 = rmse_r2.RMSECallback(solver_obj2,'../../../outputs/0.validation/discrete-4cores', E_area_centre_point1, E_area_circle1)
solver_obj2.add_callback(cb5,'timestep')


# Operation of tidal turbine farm about each turbine output through a callback
cb6 = rmse_r3.RMSECallback(solver_obj2,'../../../outputs/0.validation/discrete-4cores', E_area_centre_point2, E_area_circle2)
solver_obj2.add_callback(cb6,'timestep')

# start computer forward model
solver_obj2.iterate(update_forcings=update_forcings)
###set up interest functional and control###
power_output2= sum(cb2.average_power)
effect_two2 = (cb6.RMSEaverage + cb6.RMSEaverage)/2

power_output = (power_output1+power_output2)/2
effect_two = (effect_two1+effect_two2)/2

maxoutput, diff_uv = 3516.1413136621304, (3287.5225054565767+4479.962038685619)/2

interest_functional = (1-BE/10)*power_output/maxoutput*diff_uv- effect_two * BE/10

print(power_output1,effect_two1)
print(power_output2,effect_two2)


# specifies the control we want to vary in the optimisation
c =[Control(x) for xy in farm_options.turbine_coordinates for x in xy] + [Control(x) for x in farm_options.turbine_axis]
# print(power_output, cablecost, interest_functional)
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
    # turbines.TurbineOptimisationCallback(solver_obj, cb),
    # rmse_r2.RMSEOptimisationCallback(solver_obj,cb3),
    # turbines.TurbineOptimisationCallback(solver_obj2, cb2),
    # rmse_r3.RMSEOptimisationCallback(solver_obj2,cb6),

    # turbines.EachTurbineOptimisationCallback(solver_obj,cb2),
])
# callbacks to indicate start of forward and adjoint runs in log
def eval_cb_pre(controls):
    print_output("FORWARD RUN:")
    print_output("angle: {}".format([float(c) for c in controls]))

def derivative_cb_pre(controls):
    print_output("ADJOINT RUN:")
    print_output("angle: {}".format([float(c) for c in controls]))

# this reduces the functional J(u, m) to a function purely of the control m:
# rf(m) = J(u(m), m) where the velocities u(m) of the entire simulation
# are computed by replaying the forward model for any provided turbine coordinates m
rf = ReducedFunctional(-interest_functional, c, derivative_cb_post=callback_list,
        eval_cb_pre=eval_cb_pre, derivative_cb_pre=derivative_cb_pre)

if  0:
    # whenever the forward model is changed - for example different terms in the equation,
    # different types of boundary conditions, etc. - it is a good idea to test whether the
    # gradient computed by the adjoint is still correct, as some steps in the model may
    # not have been annotated correctly. This can be done via the Taylor test.
    # Using the standard Taylor series, we should have (for a sufficiently smooth problem):
    #   rf(td0+h*dtd) - rf(td0) - < drf/dtd(rf0), h dtd> = O(h^2)

    # we choose a random point in the control space, i.e. a randomized turbine density with
    # values between 0 and 1 and choose a random direction dtd to vary it in

    # this tests whether the above Taylor series residual indeed converges to zero at 2nd order in h as h->0

    m0 =  [Constant(x) for xy in farm_options.turbine_coordinates for x in xy] + [Constant(i) for i in farm_options.turbine_axis]
    h0 =  [Constant(1) for xy in farm_options.turbine_coordinates for x in xy] + [Constant(1) for i in farm_options.turbine_axis]

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

    lb = list(np.array([[xmin+d, ymin+d] for _ in farm_options.turbine_coordinates]).flatten()) + [0]*int(len(farm_options.turbine_axis)/2) + [180]*int(len(farm_options.turbine_axis)/2) 
    ub = list(np.array([[xmax-d, ymax-d] for _ in farm_options.turbine_coordinates]).flatten()) + [180]*int(len(farm_options.turbine_axis)/2) + [360]*int(len(farm_options.turbine_axis)/2)
    lb = [Constant(i) for i in lb]
    ub = [Constant(i) for i in ub]

    # lb = list(np.array([[xmin+d, ymin+d] for _ in farm_options.turbine_coordinates]).flatten())
    # ub = list(np.array([[xmax-d, ymax-d] for _ in farm_options.turbine_coordinates]).flatten())
    # lb = [Constant(i) for i in lb]
    # ub = [Constant(i) for i in ub]

    mdc= turbines.MinimumDistanceConstraints(farm_options.turbine_coordinates, farm_options.turbine_axis, farm_options.individual_thrust_coefficient, 40. ,optimise_layout_only)
    
    td_opt = minimize(rf, method='SLSQP', bounds=[lb,ub], constraints=mdc,
            options={'maxiter': 200, 'pgtol': 1e-3})

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
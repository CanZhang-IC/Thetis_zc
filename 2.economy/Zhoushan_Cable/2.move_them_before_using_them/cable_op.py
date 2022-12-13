from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
op2.init(log_level=INFO)
import sys
sys.path.append('../..')
import prepare.utm, prepare.myboundary
import prepare_cable.Hybrid_Code
from prepare_cable.cable_overloaded import cablelength
import numpy
import time
import yagmail
import h5py

start_time = time.time()

get_index = os.path.basename(sys.argv[0])
BE = float(get_index[:-3])
# BE = 5.0

start_from_initial = True
output_dir = '../../../outputs/2.economy/discrete/intermediate-step1-fromsame/yaw-cable-BE-'+str(BE)[:-2]
print_output(output_dir[17:])

file_dir = '../../'
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
farm_options.upwind_correction = True


xmin,ymin,xmax,ymax = 443340, 3322634, 443592, 3322848 

if start_from_initial:
    turbine_location = []
    x_space = 60
    for x in range(xmin+20,xmax-20,x_space):
        for y in range(ymin+20,ymax-20,x_space*2):
            turbine_location.append([x,y])
    for x in range(xmin+20+int(x_space/2),xmax-20,x_space):
        for y in range(ymin+20+x_space,ymax-20,x_space*2):
            turbine_location.append([x,y])
    farm_options.turbine_coordinates =[[Constant(xy[0]),Constant(xy[1])] for xy in turbine_location]

    # when 'farm_options.considering_yaw' is False, 
    # the farm_options.turbine_axis would be an empty list, 
    # so only the coordinates are optimised.
    farm_options.considering_yaw =  True
    farm_options.turbine_axis = [Constant(90) for i in range(len(farm_options.turbine_coordinates))] + [Constant(270) for i in range(len(farm_options.turbine_coordinates))]
else:
    result_output_dir = '../../../outputs/2.economy/discrete/intermediate/yaw-cable-BE-'+str(float(BE+1))[:-2]
    print('初始位置提取于',result_output_dir[17:])
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
    # farm_options.turbine_axis = [Constant(100) for i in range(len(farm_options.turbine_coordinates))] + [Constant(280) for i in range(len(farm_options.turbine_coordinates))]
    flood_dir,ebb_dir = all_controls[24:36], all_controls[36:]
    farm_options.turbine_axis = [Constant(i) for i in flood_dir] + [Constant(i) for i in ebb_dir]

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
solver_obj.load_state(492, outputdir='../../../outputs/0.validation/discrete-4cores')
#solver_obj.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=Constant(0.0))

# Operation of tidal turbine farm through a callback
cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')

# cb2 = turbines.EachTurbineFunctionalCallback(solver_obj)
# solver_obj.add_callback(cb2,'timestep')

# start computer forward model
solver_obj.iterate(update_forcings=update_forcings)

#File('turbinedensity.pvd').write(solver_obj.fields.turbine_density_2d)
###set up interest functional and control###
power_output= sum(cb.average_power)


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

c =[Control(x) for xy in farm_options.turbine_coordinates for x in xy] + [Control(i) for i in farm_options.turbine_axis]
turbine_density = Function(solver_obj.function_spaces.P1_2d, name='turbine_density')
turbine_density.interpolate(solver_obj.tidal_farms[0].turbine_density)


maxoutput, maxcost = 1754,2144

interest_functional = BE/10*power_output/maxoutput-(1-BE/10)*cablecost/maxcost
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
callback_list = optimisation.OptimisationCallbackList([
    optimisation.ConstantControlOptimisationCallback(solver_obj, array_dim=len(c)),
    optimisation.DerivativeConstantControlOptimisationCallback(solver_obj, array_dim=len(c)),
    optimisation.UserExportOptimisationCallback(solver_obj, [turbine_density, solver_obj.fields.uv_2d]),
    optimisation.FunctionalOptimisationCallback(solver_obj),
    turbines.TurbineOptimisationCallback(solver_obj, cb),
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
            options={'maxiter': 500, 'pgtol': 1e-3})

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
from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
op2.init(log_level=INFO)
from os.path import join
import detectors
import tidal_forcing
import utm
import numpy
import Hybrid_Code
from cable_overloaded import cablelength


ouput_dir = 'outputs/opcable_cost'

mesh2d = Mesh('./mesh/mesh.msh')
#timestepping options
dt = 30*60 # reduce this if solver does not converge
t_export = 30*60 
#t_end = 530*30*60 
t_end = 484*30*60 

P1 = FunctionSpace(mesh2d, "CG", 1)

# # read bathymetry code
chk = DumbCheckpoint('bathymetry', mode=FILE_READ)
bathymetry2d = Function(P1)
chk.load(bathymetry2d, name='bathymetry')
chk.close()

#read viscosity / manning boundaries code
chk = DumbCheckpoint('viscosity', mode=FILE_READ)
h_viscosity = Function(P1, name='viscosity')
chk.load(h_viscosity)
chk.close()

#manning = Function(P1,name='manning').assign(0.02)
chk = DumbCheckpoint('manning', mode=FILE_READ)
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
    x_0, y_0, utm_zone, zone_letter = utm.from_latlon(lat, 0)
    coriolis_2d = Function(FunctionSpace(mesh, 'CG', 1), name="coriolis_2d")
    coriolis_2d.interpolate(f0 + beta * (x[1] - y_0))
    return coriolis_2d
coriolis_2d =coriolis(mesh2d, 29)

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
options.element_family = "dg-cg"
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
solver_obj.bnd_functions['shallow_water'] = {
        200: {'elev': tidal_elev},  #set open boundaries to tidal_elev function
  }

def update_forcings(t):
  with timed_stage('update forcings'):
    print_output("Updating tidal field at t={}".format(t))
    elev = tidal_forcing.set_tidal_field(Function(bathymetry2d.function_space()), t)
    tidal_elev.project(elev) 
    print_output("Done updating tidal field")




# Initialise Discrete turbine farm characteristics
farm_options = DiscreteTidalTurbineFarmOptions()
farm_options.turbine_type = 'constant'
farm_options.turbine_options.thrust_coefficient = 0.6
farm_options.turbine_options.diameter = 20
farm_options.upwind_correction = False


# farm_options.turbine_coordinates =[[Constant(positions[i]), Constant(positions[i+1])] for i in range(0,len(positions),2)]
farm_options.turbine_coordinates = [[Constant(x), Constant(y)] for x in numpy.arange(443032+20, 443288-20, 60) for y in numpy.arange(3322891+20, 3323091-20, 40)]
# when 'farm_options.considering_yaw' is False, 
# the farm_options.turbine_axis would be an empty list, 
# so only the coordinates are optimised.
farm_options.considering_yaw = False
farm_options.turbine_axis = [Constant(90) for i in range(len(farm_options.turbine_coordinates))]
#add turbines to SW_equations
options.discrete_tidal_turbine_farms[2] = farm_options

# run as normal (this run will be annotated by firedrake_adjoint)
#solver_obj.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=Constant(0.0))
solver_obj.load_state(483,outputdir='./outputs/redata')

# Operation of tidal turbine farm through a callback
cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')


# start computer forward model
solver_obj.iterate(update_forcings=update_forcings)

#File('turbinedensity.pvd').write(solver_obj.fields.turbine_density_2d)
###set up interest functional and control###
power_output= sum(cb.integrated_power)

###cable length###
turbine_locations = [float(x) for xy in farm_options.turbine_coordinates for x in xy]
landpointlocation = [442500,3322750]
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank == 0:
    cableclass = Hybrid_Code.CableCostGA(turbine_locations, landpointlocation)
    order_w = cableclass.compute_cable_cost_order()
else:
    order_w = []
order_w = comm.bcast(order_w, root=0)


landpointlocation_con = [Constant(x) for x in landpointlocation]
order_con = [Constant(i) for j in order_w for i in j]
cablecost = cablelength([x for xy in farm_options.turbine_coordinates for x in xy],landpointlocation_con,order_con)
# specifies the control we want to vary in the optimisation
# optimise_angle_only = False
# if optimise_angle_only:
#     if farm_options.considering_yaw:
#         c = [Control(x) for x in farm_options.turbine_axis]
#     else:
#         raise Exception('You should turn on the yaw considering!')      
# else:
#     c = [Control(x) for xy in farm_options.turbine_coordinates for x in xy] + [Control(x) for x in farm_options.turbine_axis]

c = [Control(x) for xy in farm_options.turbine_coordinates for x in xy] 
turbine_density = Function(solver_obj.function_spaces.P1_2d, name='turbine_density')
turbine_density.interpolate(solver_obj.tidal_farms[0].turbine_density)
File('turbinedensity2.pvd').write(turbine_density)




interest_functional =power_output - cablecost
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

if 1:
    # whenever the forward model is changed - for example different terms in the equation,
    # different types of boundary conditions, etc. - it is a good idea to test whether the
    # gradient computed by the adjoint is still correct, as some steps in the model may
    # not have been annotated correctly. This can be done via the Taylor test.
    # Using the standard Taylor series, we should have (for a sufficiently smooth problem):
    #   rf(td0+h*dtd) - rf(td0) - < drf/dtd(rf0), h dtd> = O(h^2)

    # we choose a random point in the control space, i.e. a randomized turbine density with
    # values between 0 and 1 and choose a random direction dtd to vary it in

    # this tests whether the above Taylor series residual indeed converges to zero at 2nd order in h as h->0
    m1 = [[Constant(x), Constant(y)] for x in numpy.arange(443032+20, 443288-20, 60) for y in numpy.arange(3322891+20, 3323091-20, 40)]
    m0 = [i for j in m1 for i in j]

    h0 = [Constant(1) for i in range(len(farm_options.turbine_coordinates)*2)]

    minconv = taylor_test(rf, m0, h0)
    print_output("Order of convergence with taylor test (should be 2) = {}".format(minconv))

    assert minconv > 1.95

if 0:
    # Optimise the control for minimal functional (i.e. maximum profit)
    # with a gradient based optimisation algorithm using the reduced functional
    # to replay the model, and computing its derivative via the adjoint
    # By default scipy's implementation of L-BFGS-B is used, see
    #   https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html
    # options, such as maxiter and pgtol can be passed on.
    if optimise_angle_only:
        lb = [0]*len(farm_options.turbine_coordinates)
        ub = [360]*len(farm_options.turbine_coordinates)
        td_opt = minimize(rf, method='SLSQP', bounds=[lb,ub],options={'maxiter': 100, 'ptol': 1e-3})
    else:
        site_x1 = 443032.
        site_x2 = 443288.
        site_y1 = 3322891.
        site_y2 = 3323091.
        r = farm_options.turbine_options.diameter/2.

        lb = np.array([[site_x1+r, site_y1+r] for _ in farm_options.turbine_coordinates]).flatten()
        ub = np.array([[site_x2-r, site_y2-r] for _ in farm_options.turbine_coordinates]).flatten()
        
        if farm_options.considering_yaw:
            lb = list(lb) + [0]*len(farm_options.turbine_coordinates)
            ub = list(ub) + [360]*len(farm_options.turbine_coordinates)

        mdc= turbines.MinimumDistanceConstraints(farm_options.turbine_coordinates, farm_options.turbine_axis, 40.)
        
        td_opt = minimize(rf, method='SLSQP', bounds=[lb,ub], constraints=mdc,
                options={'maxiter': 5, 'pgtol': 1e-3})
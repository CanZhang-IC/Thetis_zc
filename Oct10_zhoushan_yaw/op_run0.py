from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
import numpy
op2.init(log_level=INFO)
import sys
sys.path.append('..')
import prepare.utm, prepare.myboundary_30min
# import h5py
import time

t_start = time.time()

turbine_i = 0
output_dir = '../../outputs/middle_30min-'+str(turbine_i)

mesh2d = Mesh('../mesh/mesh.msh')
#timestepping options
dt = 30*60 # reduce this if solver does not converge
t_export = 30*60 
#t_end = 1555200
#t_end = 1216800+ 13*60*60 # spring
t_end = 885600 + 13*60*60 # middle
#t_end = 612000 + 13*60*60 # neap
#t_end = 30*60

do_taylor_test = False
do_optimisation = 1

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

# result_data = []
# df = h5py.File('../../outputs/middle_30min/diagnostic_controls.hdf5','r+')
# for name, data in df.items():
#     for i in data:
#         result_data.append(i)
# print(result_data)

result_data = [
    #coordinates
    443352.00212005526, 3322646.308409128, 443352.006525044, 3322799.502752303, 443370.4567628371, 3322834.992831582, 443416.5186574596, 3322834.9959139726, 443374.8678731752, 3322766.8210353497, 443393.4028281948, 3322802.2879163884, 443440.54927117965, 3322802.996771269, 443464.6564190777, 3322834.998411545, 443417.2545525552, 3322770.4474302097, 443462.86742582766, 3322769.75287414, 443491.6995552, 3322797.5118418285, 443505.6663923886, 3322834.9991843486, 443530.14065815165, 3322725.349804717, 443545.1703302105, 3322762.4326752005, 443555.93797414476, 3322800.967399174, 443576.9362127845, 3322834.996947221,
    #axis
    99.59111587277175, 97.7452236884505, 100.3367661071047, 102.15562685125745, 98.20284924701546, 99.72541950526949, 101.95277663314398, 104.61209667439944, 100.36181525387029, 102.29380671174125, 103.5866105384266, 106.33249440850874, 102.70327494351393, 104.51145377273231, 106.97282031956618, 109.0279970709512
    ]

result_data.pop(2*turbine_i)
result_data.pop(2*turbine_i)
result_data.pop(30+turbine_i)


farm_options.turbine_coordinates = [[Constant(result_data[2*i]), Constant(result_data[2*i+1])] for i in range(int(len(result_data)/3))]
farm_options.considering_yaw = True
farm_options.turbine_axis = [Constant(i) for i in result_data[int(len(result_data)/3*2):]]
# print(len(farm_options.turbine_coordinates),len(farm_options.turbine_axis))
#add turbines to SW_equations
options.discrete_tidal_turbine_farms[2] = farm_options

###spring:676,middle:492,neap:340###
solver_obj.load_state(492, outputdir='../../outputs/redata_30min_normaldepth')
#solver_obj.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=Constant(0.0))

# Operation of tidal turbine farm through a callback
cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')


# start computer forward model

solver_obj.iterate(update_forcings=update_forcings)

#File('turbinedensity.pvd').write(solver_obj.fields.turbine_density_2d)
###set up interest functional and control###
power_output= sum(cb.integrated_power)
interest_functional = power_output

# specifies the control we want to vary in the optimisation
optimise_angle_only = True
if optimise_angle_only:
    if farm_options.considering_yaw:
        c = [Control(x) for x in farm_options.turbine_axis]
    else:
        raise Exception('You should turn on the yaw considering!')      
else:
    c = [Control(x) for xy in farm_options.turbine_coordinates for x in xy] + [Control(x) for x in farm_options.turbine_axis]

turbine_density = Function(solver_obj.function_spaces.P1_2d, name='turbine_density')
turbine_density.interpolate(solver_obj.tidal_farms[0].turbine_density)
# File('turbinedensity2.pvd').write(turbine_density)
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


if do_taylor_test:
    # whenever the forward model is changed - for example different terms in the equation,
    # different types of boundary conditions, etc. - it is a good idea to test whether the
    # gradient computed by the adjoint is still correct, as some steps in the model may
    # not have been annotated correctly. This can be done via the Taylor test.
    # Using the standard Taylor series, we should have (for a sufficiently smooth problem):
    #   rf(td0+h*dtd) - rf(td0) - < drf/dtd(rf0), h dtd> = O(h^2)

    # we choose a random point in the control space, i.e. a randomized turbine density with
    # values between 0 and 1 and choose a random direction dtd to vary it in

    # this tests whether the above Taylor series residual indeed converges to zero at 2nd order in h as h->0
    if optimise_angle_only:
        m0 = [Constant(90) for i in range(len(farm_options.turbine_coordinates))]
        h0 = [Constant(1) for i in range(len(farm_options.turbine_coordinates))]
        minconv = taylor_test(rf, m0, h0)
        print_output("Order of convergence with taylor test (should be 2) = {}".format(minconv))
        assert minconv > 1.95
    else:
        m1 = [[Constant(x), Constant(y)] for x in numpy.arange(site_x1+20, site_x2-20, 60) for y in numpy.arange(site_y1+20, site_y2-20, 40)]
        m0 = [i for j in m1 for i in j]+[Constant(90) for i in range(len(farm_options.turbine_coordinates))]
        h0 = [Constant(1) for i in range(len(farm_options.turbine_coordinates)*2)]+[Constant(1) for i in range(len(farm_options.turbine_coordinates))]

        minconv = taylor_test(rf, m0, h0)
        print_output("Order of convergence with taylor test (should be 2) = {}".format(minconv))

        assert minconv > 1.95

if do_optimisation:
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
        r = farm_options.turbine_options.diameter/2.

        lb = np.array([[site_x1+r, site_y1+r] for _ in farm_options.turbine_coordinates]).flatten()
        ub = np.array([[site_x2-r, site_y2-r] for _ in farm_options.turbine_coordinates]).flatten()
        
        if farm_options.considering_yaw:
            lb = list(lb) + [0]*len(farm_options.turbine_coordinates)
            ub = list(ub) + [360]*len(farm_options.turbine_coordinates)

        mdc= turbines.MinimumDistanceConstraints(farm_options.turbine_coordinates, farm_options.turbine_axis, 40.)
        
        td_opt = minimize(rf, method='SLSQP', bounds=[lb,ub], constraints=mdc,
                options={'maxiter': 5, 'pgtol': 1e-3})

t_end = time.time()
print('time cost:', t_end - t_start, 's')
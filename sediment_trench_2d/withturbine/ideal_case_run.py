"""
An ideal case for yaw angle optimisation
"""
# to enable a gradient-based optimisation using the adjoint to compute
# gradients, we need to import from thetis_adjoint instead of thetis. This
# ensure all firedrake operations in the Thetis model are annotated
# automatically, in such a way that we can rerun the model with different input
# parameters, and also derive the adjoint-based gradient of a specified input
# (the functional) with respect to a specified input (the control)
from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
op2.init(log_level=INFO)
import sys
sys.path.append('..')
import os
import numpy
import time
import h5py
import yagmail
import matplotlib.pyplot as plt
import bathymetry_changes_callback 

t_start = time.time()

conservative = False

P_factor = 0.5
H = 40
distance = 10
speed = 2
output_dir = '../../../outputs/sediment/test/run_100_12_op_0.5'

### set up the Thetis solver obj as usual ##
mesh2d = Mesh('../../prepare_ideal_meshes/conference_mesh2_with_effected_area.msh')

morfac = 100
dt = 30
end_time = 30*1200

diffusivity = 0.15

#set viscosity bumps at in-flow boundaries.
P1_2d = FunctionSpace(mesh2d, 'CG', 1)
V = FunctionSpace(mesh2d, "CG", 1)
DG_2d = FunctionSpace(mesh2d, "DG", 1)
vector_dg = VectorFunctionSpace(mesh2d, "DG", 1)

h_viscosity = Constant(1e-6)

bathymetry_2d = Function(V, name='bathymetry_2d').assign(Constant(40))
initial_bathymetry_2d = Function(V, name='bathymetry_2d').assign(Constant(40))

# initialise velocity and elevation
chk = DumbCheckpoint("../../../outputs/sediment/hydrodynamics_trench/elevation", mode=FILE_READ)
elev = Function(DG_2d, name="elevation")
chk.load(elev)
chk.close()

chk = DumbCheckpoint('../../../outputs/sediment/hydrodynamics_trench/velocity', mode=FILE_READ)
uv = Function(vector_dg, name="velocity")
chk.load(uv)
chk.close()

# create solver and set options
solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry_2d)
options = solver_obj.options

options.sediment_model_options.solve_suspended_sediment = True
options.sediment_model_options.use_bedload = True
options.sediment_model_options.solve_exner = True

options.sediment_model_options.use_sediment_conservative_form = conservative
options.sediment_model_options.average_sediment_size = Constant(160*(10**(-6)))
options.sediment_model_options.bed_reference_height = Constant(0.025)
options.sediment_model_options.morphological_acceleration_factor = Constant(morfac)

options.simulation_end_time = end_time/morfac
options.simulation_export_time = options.simulation_end_time/45
options.output_directory = output_dir
options.check_volume_conservation_2d = True

if options.sediment_model_options.solve_suspended_sediment:
    options.fields_to_export = ['sediment_2d', 'uv_2d', 'elev_2d', 'bathymetry_2d']  # note exporting bathymetry must be done through export func
    options.sediment_model_options.check_sediment_conservation = True
else:
    options.fields_to_export = ['uv_2d', 'elev_2d', 'bathymetry_2d']  # note exporting bathymetry must be done through export func

options.element_family = 'dg-dg'
options.timestepper_type = 'CrankNicolson'
options.timestepper_options.implicitness_theta = 1.0
options.timestepper_options.use_semi_implicit_linearization = True
# using direct solver as PressurePicard does not work with dolfin-adjoint (due to .split() not being annotated correctly)
options.timestepper_options.solver_parameters = {'snes_monitor': None,
                                                'snes_rtol': 1e-5,
                                                'ksp_type': 'preonly',
                                                'pc_type': 'lu',
                                                'pc_factor_mat_solver_type': 'mumps',
                                                'mat_type': 'aij'
                                                }
options.norm_smoother = Constant(0.1)

# using nikuradse friction
options.nikuradse_bed_roughness = Constant(3*options.sediment_model_options.average_sediment_size)
options.horizontal_viscosity = h_viscosity
options.horizontal_diffusivity = Constant(diffusivity)

if not hasattr(options.timestepper_options, 'use_automatic_timestep'):
    options.timestep = dt

# assign boundary conditions
left_tag = 1
right_tag = 2
coasts_tag = 3
tidal_elev = Function(get_functionspace(mesh2d, "CG", 1), name='tidal_elev').assign(0.0)

tidal_vel = Function(V).assign(0.0)

solver_obj.bnd_functions['shallow_water'] = {
    right_tag: {'uv':as_vector((-Constant(speed), 0.0))},
    left_tag: {'elev': tidal_elev},
}

# Initialise Discrete turbine farm characteristics
farm_options = DiscreteTidalTurbineFarmOptions()
farm_options.turbine_type = 'constant'
farm_options.turbine_options.thrust_coefficient = 0.6
farm_options.turbine_options.diameter = 20
farm_options.upwind_correction = True

turbine_location = []
# for i in range(850,1200,200):
#     for j in range(250, 400, 50):
#         turbine_location.append([i,j])
# for i in range(950,1200,200):
#     for j in range(275, 400, 50):
#         turbine_location.append([i,j])

for i in range(1450,1800,100):
    for j in range(250, 400, 50):
        turbine_location.append([i,j])

# turbine_location = [1533.7677940147619, 238.42127956411298, 1523.1986992983088, 304.22227845195505, 1508.7693318280867, 341.5290252202772, 1569.6309058898694, 256.13673063855407, 1558.8645630099643, 322.33156242300277, 1578.8926348489072, 356.9563561204352, 1644.677449632057, 224.80173130452548, 1592.4992077269844, 288.9550299718026, 1615.3308858611126, 373.45647319532765, 1751.5742443339452, 210.0, 1738.7286845818496, 273.5534701577489, 1651.7494487242504, 390.0]

farm_options.turbine_coordinates =[
    [Constant(xy[0]),Constant(xy[1])] for xy in turbine_location
    # [Constant(turbine_location[2*i]),Constant(turbine_location[2*i+1])] for i in range(int(len(turbine_location)/2))
    ]
farm_options.considering_yaw = False

#add turbines to SW_equations
options.discrete_tidal_turbine_farms[2] = farm_options

# def update_forcings(t):
#     print_output("Updating tidal elevation at t = {}".format(t))
#     tidal_vel.project(-Constant(speed))

#set initial condition
solver_obj.assign_initial_conditions(uv=as_vector((-Constant(speed), 0.0)), elev=tidal_elev)

# Operation of tidal turbine farm through a callback
cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')

#Effected area location
E_area_centre_point = [1400,250]
E_area_circle = 40

# Operation of tidal turbine farm about each turbine output through a callback
cb2 = bathymetry_changes_callback.BathymetryCallback(solver_obj, initial_bathymetry_2d, E_area_centre_point, E_area_circle)
solver_obj.add_callback(cb2,'timestep')

# start computer forward model
solver_obj.iterate()#update_forcings=update_forcings)

# ###set up interest functional and control###
power_output= sum(cb.average_power)
maxoutput, maxeffect = 7579.476290597467, 994.0157143521978
interest_functional = (P_factor*(power_output/maxoutput)-(1-P_factor)*(cb2.b_c_average/maxeffect))*maxoutput
print(interest_functional,power_output,cb2.b_c_average)


# specifies the control we want to vary in the optimisation
optimise_angle_only = False
optimise_layout_only = True
# if optimise_angle_only:
#     if farm_options.considering_yaw:
#         c = [Control(x) for x in farm_options.turbine_axis]
#     else:
#         raise Exception('You should turn on the yaw considering!')    
# elif optimise_layout_only:

#     c =  [Control(x) for xy in farm_options.turbine_coordinates for x in xy] 
# else:
#     c = [Control(x) for xy in farm_options.turbine_coordinates for x in xy] + [Control(x) for x in farm_options.turbine_axis]
c = [Control(x) for xy in farm_options.turbine_coordinates for x in xy] #+  [Control(x) for x in farm_options.turbine_axis] + [Control(x) for x in farm_options.individual_thrust_coefficient]

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
    # m0 =  [Constant(x) for xy in farm_options.turbine_coordinates for x in xy] + [Constant(x) for x in farm_options.turbine_axis]
    # h0 =  [Constant(1) for xy in farm_options.turbine_coordinates for x in xy] +[Constant(1) for x in farm_options.turbine_axis]
    m0 =  [Constant(x) for xy in farm_options.turbine_coordinates for x in xy]#+[Constant(x) for x in farm_options.turbine_axis] + [Constant(x) for x in farm_options.individual_thrust_coefficient] 
    h0 =  [Constant(1) for xy in farm_options.turbine_coordinates for x in xy] #+ [Constant(1) for x in farm_options.turbine_axis] + [Constant(0.1) for x in farm_options.individual_thrust_coefficient] 
    # m0 =  [Constant(x) for xy in farm_options.turbine_coordinates for x in xy] 
    # h0 =  [Constant(1) for xy in farm_options.turbine_coordinates for x in xy] 
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
    if optimise_angle_only:
        lb = [-90]*len(farm_options.turbine_axis) + [0.1]*len(farm_options.individual_thrust_coefficient)
        ub = [90]*len(farm_options.turbine_axis) + [0.6]*len(farm_options.individual_thrust_coefficient)
        td_opt = minimize(rf, method='SLSQP', bounds=[lb,ub],options={'maxiter': 200, 'ptol': 1e-3})
    elif optimise_layout_only:
        site_x = 400.
        site_y = 200.
        site_x_start = 1400.
        site_y_start = 200.
        r = farm_options.turbine_options.diameter/2.

        lb = np.array([[site_x_start+r, site_y_start+r] for _ in farm_options.turbine_coordinates]).flatten()
        ub = np.array([[site_x_start+site_x-r, site_y_start+site_y-r] for _ in farm_options.turbine_coordinates]).flatten()

        mdc= turbines.MinimumDistanceConstraints(farm_options.turbine_coordinates, farm_options.turbine_axis, farm_options.individual_thrust_coefficient, 40. ,optimise_layout_only)
        
        td_opt = minimize(rf, method='SLSQP', bounds=[lb,ub], constraints=mdc,
                options={'maxiter': 200, 'pgtol': 1e-3})
    else:
        site_x = 400.
        site_y = 200.
        site_x_start = 1400.
        site_y_start = 200.
        r = farm_options.turbine_options.diameter/2.

        lb = np.array([[site_x_start+r, site_y_start+r] for _ in farm_options.turbine_coordinates]).flatten()
        ub = np.array([[site_x_start+site_x-r, site_y_start+site_y-r] for _ in farm_options.turbine_coordinates]).flatten()
        
        if farm_options.considering_yaw:
            lb = list(lb) + [-90]*len(farm_options.turbine_axis)
            ub = list(ub) + [90]*len(farm_options.turbine_axis)

        if farm_options.considering_individual_thrust_coefficient:
            lb = list(lb) + [0.1]*len(farm_options.individual_thrust_coefficient)
            ub = list(ub) + [0.6]*len(farm_options.individual_thrust_coefficient)

        mdc= turbines.MinimumDistanceConstraints(farm_options.turbine_coordinates, farm_options.turbine_axis, farm_options.individual_thrust_coefficient, 40. ,optimise_layout_only)
        
        td_opt = minimize(rf, method='SLSQP', bounds=[lb,ub], constraints=mdc,
                options={'maxiter': 200, 'pgtol': 1e-3})

t_end = time.time()
print('time cost: {0:.2f}min'.format((t_end - t_start)/60))

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank == 0 :
    yag = yagmail.SMTP(user = '623001493@qq.com',password = 'Zc623oo1493', host = 'smtp.qq.com')
    yag.send(to = ['canzhang2019@gmail.com'], subject = 'Python done', contents = ['College computer Python Finished'+str(P_factor)])
else:
    pass


# final_bathymetry = Function(V, name = 'final_bathymetry').interpolate(-solver_obj.fields.bathymetry_2d) 
# File('final_bathymetry.pvd').write(final_bathymetry)

# bathymetry_changes = Function(V, name='bathymetry_changes').interpolate(initial_bathymetry_2d + final_bathymetry)
# File('bathymetry_changes.pvd').write(bathymetry_changes)

# # record final bathymetry for plotting
# xaxisthetis1 = []
# baththetis1 = []

# for i in np.linspace(100, 900, 1000):
#     xaxisthetis1.append(i)
#     if conservative:
#         baththetis1.append(-solver_obj.fields.bathymetry_2d.at([i, 50]))
#     else:
#         baththetis1.append(-solver_obj.fields.bathymetry_2d.at([i, 50]))

# if os.getenv('THETIS_REGRESSION_TEST') is None:
#     # Compare model and experimental results
#     # (this part is skipped when run as a test)
#     # data = np.genfromtxt('experimental_data.csv', delimiter=',')

#     # plt.scatter([i[0] for i in data], [i[1] for i in data], label='Experimental Data')

    # plt.plot(xaxisthetis1, baththetis1, label='Thetis')
    # plt.legend()
    # plt.show()




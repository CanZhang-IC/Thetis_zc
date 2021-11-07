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
import rmse_r2
import numpy
import time

t_start = time.time()

# get_index = os.path.basename(sys.argv[0])
# namelength = len('closed_boundary_upwind')
# P_factor = float(get_index[namelength:-3])

# if angle_H[0] == '0':
#     angle = 0
#     H = int(angle_H[1:])
# else:
#     angle = int(angle_H[:2])
#     H = int(angle_H[2:])

# P_factor = 1
angle, H = 0, 40 

speed = 2
output_dir = '../../../outputs/6.yaw_environment/Yaw_Ideal/restart-conference_mesh2_with_effected_area-one-core'
#Comment for testing forward model

### set up the Thetis solver obj as usual ##
mesh2d = Mesh('../../prepare_ideal_meshes/conference_mesh2_with_effected_area.msh')

timestep = 60
t_export = 2 * timestep
t_end = 50*t_export #12000


#set viscosity bumps at in-flow boundaries.
P1_2d = FunctionSpace(mesh2d, 'CG', 1)
x = SpatialCoordinate(mesh2d)
h_viscosity = Constant(1e-3)
# h_viscosity = Function(P1_2d).interpolate(conditional(le(x[0], 20), 20+viscosity-x[0], conditional(ge(x[0],1980),x[0]-viscosity,viscosity)))
# File(output_dir+'/viscosity.pvd').write(h_viscosity)

# create solver and set options
solver_obj = solver2d.FlowSolver2d(mesh2d, Constant(H))
options = solver_obj.options
options.timestep = timestep
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.output_directory = output_dir
options.check_volume_conservation_2d = True
options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']
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
options.horizontal_viscosity = h_viscosity
options.quadratic_drag_coefficient = Constant(0.0025)

# assign boundary conditions
left_tag = 1
right_tag = 2
coasts_tag = 3
tidal_elev = Function(get_functionspace(mesh2d, "CG", 1), name='tidal_elev')

tidal_vel = Function(P1_2d).assign(0.0)

solver_obj.bnd_functions['shallow_water'] = {
    right_tag: {'un':tidal_vel},
    left_tag: {'elev': tidal_elev},
}

# # a function to update the tidal_elev bc value every timestep
# x = SpatialCoordinate(mesh2d)

# # Initialise Discrete turbine farm characteristics
# farm_options = DiscreteTidalTurbineFarmOptions()
# farm_options.turbine_type = 'constant'
# farm_options.turbine_options.thrust_coefficient = 0.6
# farm_options.turbine_options.diameter = 20
# farm_options.turbine_options.C_support = 1
# farm_options.turbine_options.A_support = H/2
# farm_options.upwind_correction = True

# turbine_location = []

# for i in range(850,1200,200):
#     for j in range(250, 400, 50):
#         turbine_location.append([i,j])
# for i in range(950,1200,200):
#     for j in range(275, 400, 50):
#         turbine_location.append([i,j])

# for i in range(1450,1800,100):
#     for j in range(250, 400, 50):
#         turbine_location.append([i,j])

# farm_options.turbine_coordinates =[
#     [Constant(xy[0]),Constant(xy[1])] for xy in turbine_location
#     # [Constant(600),Constant(250)]
#     ]

# farm_options.considering_yaw = True
# farm_options.turbine_axis = [Constant(angle) for i in range(len(farm_options.turbine_coordinates))]

# farm_options.considering_individual_thrust_coefficient = True
# farm_options.individual_thrust_coefficient = [Constant(0.6) for i in range(len(farm_options.turbine_coordinates))]

# farm_options.farm_alpha = Function(P1_2d)
#add turbines to SW_equations
# options.discrete_tidal_turbine_farms[2] = farm_options

def update_forcings(t):
    print_output("Updating tidal elevation at t = {}".format(t))
    tidal_vel.project(-Constant(speed))

    # uv, eta = split(solver_obj.fields.solution_2d)
    # alpha_flood = sum((conditional((x[0]-xi)**2+(x[1]-yi)**2 < farm_options.turbine_options.diameter**2, alphai/180*pi, 0) \
    #             for alphai,(xi,yi) in zip(farm_options.turbine_axis[:len(farm_options.turbine_coordinates)],farm_options.turbine_coordinates)))
    # alpha_ebb   = sum((conditional((x[0]-xi)**2+(x[1]-yi)**2 < farm_options.turbine_options.diameter**2, alphai/180*pi, 0) \
    #             for alphai,(xi,yi) in zip(farm_options.turbine_axis[len(farm_options.turbine_coordinates):],farm_options.turbine_coordinates)))
    # alphainterpolation = conditional(uv[0] > 0, alpha_flood, alpha_ebb)
    # farm_options.farm_alpha.interpolate(alphainterpolation)

solver_obj.assign_initial_conditions(uv=as_vector((-Constant(speed), 0.0)), elev=tidal_elev)

# # Operation of tidal turbine farm through a callback
# cb = turbines.TurbineFunctionalCallback(solver_obj)
# solver_obj.add_callback(cb, 'timestep')

# #Effected area location
# E_area_centre_point = [450,250]
# E_area_circle = 40

# # Operation of tidal turbine farm about each turbine output through a callback
# cb2 = rmse_r2.RMSECallback(solver_obj,'../../../outputs/6.yaw_environment/Yaw_Ideal/test_yaw_effect_restart', E_area_centre_point, E_area_circle)
# solver_obj.add_callback(cb2,'timestep')

# start computer forward model
solver_obj.iterate(update_forcings=update_forcings)

# # ###set up interest functional and control###
# power_output= sum(cb.current_power)
# interest_functional = (P_factor*(power_output/7643.52)-(1-P_factor)*(cb2.RMSE_current[-1]/7040.53))*7643.52
# print(interest_functional,power_output,cb2.RMSE_current[-1])


# # specifies the control we want to vary in the optimisation
# optimise_angle_only = True
# optimise_layout_only = False
# # if optimise_angle_only:
# #     if farm_options.considering_yaw:
# #         c = [Control(x) for x in farm_options.turbine_axis]
# #     else:
# #         raise Exception('You should turn on the yaw considering!')    
# # elif optimise_layout_only:

# #     c =  [Control(x) for xy in farm_options.turbine_coordinates for x in xy] 
# # else:
# #     c = [Control(x) for xy in farm_options.turbine_coordinates for x in xy] + [Control(x) for x in farm_options.turbine_axis]
# c = [Control(x) for x in farm_options.individual_thrust_coefficient]
# turbine_density = Function(solver_obj.function_spaces.P1_2d, name='turbine_density')
# turbine_density.interpolate(solver_obj.tidal_farms[0].turbine_density)
# # a number of callbacks to provide output during the optimisation iterations:
# # - ControlsExportOptimisationCallback export the turbine_friction values (the control)
# #            to outputs/control_turbine_friction.pvd. This can also be used to checkpoint
# #            the optimisation by using the export_type='hdf5' option.
# # - DerivativesExportOptimisationCallback export the derivative of the functional wrt
# #            the control as computed by the adjoint to outputs/derivative_turbine_friction.pvd
# # - UserExportOptimisationCallback can be used to output any further functions used in the
# #            forward model. Note that only function states that contribute to the functional are
# #            guaranteed to be updated when the model is replayed for different control values.
# # - FunctionalOptimisationCallback simply writes out the (scaled) functional values
# # - the TurbineOptimsationCallback outputs the average power, cost and profit (for each
# #            farm if multiple are defined)
# callback_list = optimisation.OptimisationCallbackList([
#     optimisation.ConstantControlOptimisationCallback(solver_obj, array_dim=len(c)),
#     optimisation.DerivativeConstantControlOptimisationCallback(solver_obj, array_dim=len(c)),
#     optimisation.UserExportOptimisationCallback(solver_obj, [turbine_density, solver_obj.fields.uv_2d]),
#     optimisation.FunctionalOptimisationCallback(solver_obj),
# ])
# # callbacks to indicate start of forward and adjoint runs in log
# def eval_cb_pre(controls):
#     print_output("FORWARD RUN:")
#     print_output("angle: {}".format([float(c) for c in controls]))

# def derivative_cb_pre(controls):
#     print_output("ADJOINT RUN:")
#     print_output("angle: {}".format([float(c) for c in controls]))

# # this reduces the functional J(u, m) to a function purely of the control m:
# # rf(m) = J(u(m), m) where the velocities u(m) of the entire simulation
# # are computed by replaying the forward model for any provided turbine coordinates m
# rf = ReducedFunctional(-interest_functional, c, derivative_cb_post=callback_list,
#         eval_cb_pre=eval_cb_pre, derivative_cb_pre=derivative_cb_pre)

# if 0:
#     # whenever the forward model is changed - for example different terms in the equation,
#     # different types of boundary conditions, etc. - it is a good idea to test whether the
#     # gradient computed by the adjoint is still correct, as some steps in the model may
#     # not have been annotated correctly. This can be done via the Taylor test.
#     # Using the standard Taylor series, we should have (for a sufficiently smooth problem):
#     #   rf(td0+h*dtd) - rf(td0) - < drf/dtd(rf0), h dtd> = O(h^2)

#     # we choose a random point in the control space, i.e. a randomized turbine density with
#     # values between 0 and 1 and choose a random direction dtd to vary it in

#     # this tests whether the above Taylor series residual indeed converges to zero at 2nd order in h as h->0
#     # m0 =  [Constant(x) for xy in farm_options.turbine_coordinates for x in xy] + [Constant(x) for x in farm_options.turbine_axis]
#     # h0 =  [Constant(1) for xy in farm_options.turbine_coordinates for x in xy] +[Constant(1) for x in farm_options.turbine_axis]
#     m0 =  [Constant(0.6) for x in farm_options.turbine_axis]
#     h0 =  [Constant(0.1) for x in farm_options.turbine_axis]
#     # m0 =  [Constant(x) for xy in farm_options.turbine_coordinates for x in xy] 
#     # h0 =  [Constant(1) for xy in farm_options.turbine_coordinates for x in xy] 
#     minconv = taylor_test(rf, m0, h0)
#     print_output("Order of convergence with taylor test (should be 2) = {}".format(minconv))

#     assert minconv > 1.95

# if 0:
#     # Optimise the control for minimal functional (i.e. maximum profit)
#     # with a gradient based optimisation algorithm using the reduced functional
#     # to replay the model, and computing its derivative via the adjoint
#     # By default scipy's implementation of L-BFGS-B is used, see
#     #   https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html
#     # options, such as maxiter and pgtol can be passed on.
#     if optimise_angle_only:
#         lb = [-90]*len(farm_options.turbine_axis)
#         ub = [90]*len(farm_options.turbine_axis)
#         td_opt = minimize(rf, method='SLSQP', bounds=[lb,ub],options={'maxiter': 200, 'ptol': 1e-3})
#     elif optimise_layout_only:
#         site_x = 400.
#         site_y = 200.
#         site_x_start = 800.
#         site_y_start = 200.
#         r = farm_options.turbine_options.diameter/2.

#         lb = np.array([[site_x_start+r, site_y_start+r] for _ in farm_options.turbine_coordinates]).flatten()
#         ub = np.array([[site_x_start+site_x-r, site_y_start+site_y-r] for _ in farm_options.turbine_coordinates]).flatten()

#         mdc= turbines.MinimumDistanceConstraints(farm_options.turbine_coordinates, farm_options.turbine_axis, 40. ,optimise_layout_only)
        
#         td_opt = minimize(rf, method='SLSQP', bounds=[lb,ub], constraints=mdc,
#                 options={'maxiter': 100, 'pgtol': 1e-3})
#     else:
#         site_x = 400.
#         site_y = 200.
#         site_x_start = 800.
#         site_y_start = 200.
#         r = farm_options.turbine_options.diameter/2.

#         lb = np.array([[site_x_start+r, site_y_start+r] for _ in farm_options.turbine_coordinates]).flatten()
#         ub = np.array([[site_x_start+site_x-r, site_y_start+site_y-r] for _ in farm_options.turbine_coordinates]).flatten()
        
#         if farm_options.considering_yaw:
#             lb = list(lb) + [-90]*len(farm_options.turbine_axis)
#             ub = list(ub) + [90]*len(farm_options.turbine_axis)

#         mdc= turbines.MinimumDistanceConstraints(farm_options.turbine_coordinates, farm_options.turbine_axis, 40. ,optimise_layout_only)
        
#         td_opt = minimize(rf, method='SLSQP', bounds=[lb,ub], constraints=mdc,
#                 options={'maxiter': 100, 'pgtol': 1e-3})

# t_end = time.time()
# print('time cost: {0:.2f}min'.format((t_end - t_start)/60))

# # comm = MPI.COMM_WORLD
# # rank = comm.Get_rank()
# # if rank == 0 :
# #     yag = yagmail.SMTP(user = '623001493@qq.com',password = 'Zc623oo1493', host = 'smtp.qq.com')
# #     yag.send(to = ['canzhang2019@gmail.com'], subject = 'Python done', contents = ['Python Finished'])
# # else:
# #     pass

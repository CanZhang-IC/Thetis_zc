"""
An ideal case for yaw angle optimisation
"""
# to enable a gradient-based optimisation using the adjoint to compute
# gradients, we need to import from thetis_adjoint instead of thetis. This
# ensure all firedrake operations in the Thetis model are annotated
# automatically, in such a way that we can rerun the model with different input
# parameters, and also derive the adjoint-based gradient of a specified input
# (the functional) with respect to a specified input (the control)
from ast import Constant
from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
op2.init(log_level=INFO)
import sys
sys.path.append('..')
import os
import time

import yagmail
import matplotlib.pyplot as plt
import numpy as np
import h5py


t_start = time.time()

get_index = os.path.basename(sys.argv[0])
the_fac = float(get_index[:-3])

# the_fac = 16
myparameters={}
myparameters['conservative'] = False
myparameters['speed'] = 2 * the_fac
myparameters['run_time'] = 1*2
myparameters['morfac'] = 100
myparameters['bath_value'] = 40
myparameters['sediment_size'] = 16*1e-5
myparameters['bnd_sediment'] = 0
myparameters['diffusivity'] = 0.15
myparameters['h_viscosity'] = 1e-6
myparameters['mor_viscosity'] = 1e-6
myparameters['bed_reference_height'] = 3 * myparameters['sediment_size']
myparameters['porosity'] = 0.4
myparameters['bedload'] = True
myparameters['suspended'] = True

# choose directory to output results
st = datetime.datetime.fromtimestamp(t_start).strftime('%Y-%m-%d')# %H:%M:%S')
special_mark = 'inflow_velocity/'+str(float(the_fac))
output_dir = './outputs/'+special_mark

with open(output_dir+'/myparameters.txt','a+') as f:
    f.write(str(myparameters))


### set up the Thetis solver obj as usual ##
mesh2d = Mesh('./mesh.msh')

dt = 30
end_time = 3600*24*myparameters['run_time']

#set viscosity bumps at in-flow boundaries.
P1_2d = FunctionSpace(mesh2d, 'CG', 1)
V = FunctionSpace(mesh2d, "CG", 1)
DG_2d = FunctionSpace(mesh2d, "DG", 1)
vector_dg = VectorFunctionSpace(mesh2d, "DG", 1)


bathymetry_2d = Function(V, name='bathymetry_2d').assign(Constant(myparameters['bath_value']))
# chk = DumbCheckpoint('/home/can/Git_thetis/outputs/sediment/test/'+'M100-H40-24h-2022-03-03-outputs-final_bathymetry',mode=FILE_READ)
# bathymetry_2d = Function(V, name='final_bathymetry')
# chk.load(bathymetry_2d)
# chk.close()

initial_bathymetry_2d = Function(V, name='ini_bathymetry_2d').assign(Constant(myparameters['bath_value']))

# initialise velocity and elevation
chk = DumbCheckpoint("./outputs/pure-hydro-v1/elevation", mode=FILE_READ)
elev = Function(DG_2d, name="elevation")
chk.load(elev)
chk.close()

chk = DumbCheckpoint('./outputs/pure-hydro-v1/velocity', mode=FILE_READ)
uv = Function(vector_dg, name="velocity")
chk.load(uv)
chk.close()

# create solver and set options
solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry_2d)
options = solver_obj.options

options.sediment_model_options.solve_suspended_sediment = myparameters['suspended']
options.sediment_model_options.use_bedload = myparameters['bedload']
options.sediment_model_options.solve_exner = True
# options.sediment_model_options.use_angle_correction = True
# options.sediment_model_options.use_slope_mag_correction = True
# options.sediment_model_options.use_secondary_current = True
# options.sediment_model_options.use_advective_velocity_correction = False

options.sediment_model_options.use_sediment_conservative_form = myparameters['conservative']
options.sediment_model_options.average_sediment_size = Constant(myparameters['sediment_size'])
options.sediment_model_options.bed_reference_height = Constant(myparameters['bed_reference_height'])
options.sediment_model_options.morphological_acceleration_factor = Constant(myparameters['morfac'])
options.sediment_model_options.morphological_viscosity = Constant(myparameters['mor_viscosity'])
options.sediment_model_options.porosity = Constant(myparameters['porosity'])
options.horizontal_viscosity = Constant(myparameters['h_viscosity'])
options.nikuradse_bed_roughness = Constant(3*options.sediment_model_options.average_sediment_size)
options.horizontal_diffusivity = Constant(myparameters['diffusivity'])

options.simulation_end_time = end_time/myparameters['morfac']
options.simulation_export_time = dt#options.simulation_end_time/45
options.output_directory = output_dir
options.check_volume_conservation_2d = True

if options.sediment_model_options.solve_suspended_sediment:
    options.fields_to_export = ['sediment_2d', 'uv_2d', 'elev_2d', 'bathymetry_2d']  # note exporting bathymetry must be done through export func
    options.sediment_model_options.check_sediment_conservation = True
else:
    options.fields_to_export = ['uv_2d', 'elev_2d', 'bathymetry_2d']  # note exporting bathymetry must be done through export func
options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d', 'bathymetry_2d']
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


if not hasattr(options.timestepper_options, 'use_automatic_timestep'):
    options.timestep = dt

# assign boundary conditions
left_tag = 1
right_tag = 2
coasts_tag = 3

tidal_elev = Function(get_functionspace(mesh2d, "CG", 1), name='tidal_elev').assign(0.0)


solver_obj.bnd_functions['shallow_water'] = {
    right_tag: {'flux':Constant(-myparameters['speed']*40*300)},
    left_tag: {'elev': tidal_elev},
}

# Initialise Discrete turbine farm characteristics
farm_options = DiscreteTidalTurbineFarmOptions()
farm_options.turbine_type = 'constant'
farm_options.turbine_options.thrust_coefficient = 0.6
farm_options.turbine_options.diameter = 20
farm_options.upwind_correction = True

turbine_location = []
# df_file = h5py.File('./outputs/turbine-optimisation-nosediment/2022-03-08-original/diagnostic_controls.hdf5','r+')
# for name,data in df_file.items():
#     for i in data[-1]:
#         turbine_location.append(i)

# for i in range(650,800,80):
#     for j in range(80, 220, 40):
#         turbine_location.append([i,j])
# for i in range(690,800,80):
#     for j in range(100, 230, 40):
#         turbine_location.append([i,j])

# for i in range(650,800,40):
#     for j in range(80, 220, 40):
#         turbine_location.append([i,j])

farm_options.turbine_coordinates =[
    # [Constant(xy[0]),Constant(xy[1])] for xy in turbine_location
    # [Constant(turbine_location[2*i]),Constant(turbine_location[2*i+1])] for i in range(int(len(turbine_location)/2))
    [Constant(700),Constant(150)]
    ]
farm_options.considering_yaw = False

#add turbines to SW_equations
options.discrete_tidal_turbine_farms[2] = farm_options

# def update_forcings(t):
#     print_output("Updating tidal elevation at t = {}".format(t))
#     tidal_vel.project(-Constant(speed))

#set initial condition
if options.sediment_model_options.solve_suspended_sediment:
    solver_obj.bnd_functions['sediment'] = {
        right_tag:{'flux':Constant(-myparameters['bnd_sediment']*40*300),'equilibrium':None},
        left_tag:{'elev': tidal_elev}
    }
    solver_obj.assign_initial_conditions(uv=uv, elev=elev)
else:
    solver_obj.assign_initial_conditions(uv=uv, elev=elev)

# Operation of tidal turbine farm through a callback
cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')

# #Effected area location
# E_area_centre_point = [400,250]
# E_area_circle = 40

# # Operation of tidal turbine farm about each turbine output through a callback
# cb2 = bathymetry_changes_callback.BathymetryCallback(solver_obj, initial_bathymetry_2d, E_area_centre_point, E_area_circle)
# solver_obj.add_callback(cb2,'timestep')

# start computer forward model
solver_obj.iterate()#update_forcings=update_forcings)


# ###set up interest functional and control###
power_output= sum(cb.current_power)
# maxoutput, maxeffect = 7579.476290597467, 994.0157143521978
interest_functional = power_output#(P_factor*(power_output/maxoutput)-(1-P_factor)*(cb2.b_c_average/maxeffect))*maxoutput
print(power_output)



# specifies the control we want to vary in the optimisation
optimise_angle_only = False
optimise_layout_only = True

c = [Control(x) for xy in farm_options.turbine_coordinates for x in xy] 

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

if 0:
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
        site_x = 170.
        site_y = 170.
        site_x_start = 640.
        site_y_start = 70.
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

        lb = np.array([[650, 80] for _ in farm_options.turbine_coordinates]).flatten()
        ub = np.array([[800, 220] for _ in farm_options.turbine_coordinates]).flatten()
        
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

if the_fac == 50:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0 :
        yag = yagmail.SMTP(user = '623001493@qq.com',password = 'ouehigyjxpidbbcj', host = 'smtp.qq.com')
        yag.send(to = ['623001493@qq.com'], subject = 'Python done', contents = ['College computer Python Finished'])
    else:
        pass

# # final_bathymetry = Function(V, name = 'final_bathymetry').interpolate(-solver_obj.fields.bathymetry_2d) 
# # File('final_bathymetry.pvd').write(final_bathymetry)

# # bathymetry_changes = Function(V, name='bathymetry_changes').interpolate(initial_bathymetry_2d + final_bathymetry)
# # File('bathymetry_changes.pvd').write(bathymetry_changes)

# # # record final bathymetry for plotting
# # xaxisthetis1 = []
# # baththetis1 = []

# # for i in np.linspace(100, 900, 1000):
# #     xaxisthetis1.append(i)
# #     if conservative:
# #         baththetis1.append(-solver_obj.fields.bathymetry_2d.at([i, 300]))
# #     else:
# #         baththetis1.append(-solver_obj.fields.bathymetry_2d.at([i, 300]))

# # if os.getenv('THETIS_REGRESSION_TEST') is None:
# #     # Compare model and experimental results
# #     # (this part is skipped when run as a test)
# #     # data = np.genfromtxt('experimental_data.csv', delimiter=',')

# #     # plt.scatter([i[0] for i in data], [i[1] for i in data], label='Experimental Data')

# #     plt.plot(xaxisthetis1, baththetis1, label='Thetis')
# #     plt.legend()
# #     plt.savefig('./fig-'+special_mark+'.png',dpi=300)

# # with open('power.txt','a+') as f:
# #     f.write(str(morfac)+'\t')
# #     f.write(str(interest_functional)+'\t')
# #     f.write(str(run_time)+'\n')



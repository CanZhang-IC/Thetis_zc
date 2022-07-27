from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
import numpy
op2.init(log_level=INFO)
import sys
sys.path.append('../..')
import prepare_continuous.utm, prepare_continuous.myboundary, prepare_continuous.detectors
import os
import time
import yagmail

start_time = time.time()

get_index = os.path.basename(sys.argv[0])
BE = float(get_index[:-3])

output_dir = '../../../outputs/2.economy/continuous/intermediate/BE'+str(BE)[:-2]
print_output(output_dir[17:])

file_dir = '../../'
mesh2d = Mesh(file_dir+'mesh_continuous/mesh.msh')
#timestepping options
dt = 5*60 # reduce this if solver does not converge
t_export = 30*60 
# t_end = 1555200
# t_end = 1216800+ 13*60*60 # spring
t_end = 885600 + 13*60*60 # middle
# t_end = 612000 + 13*60*60 # neap

test_gradient, optimise = 0,1

turbine_area_PhyID = 2

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

#manning = Function(P1,name='manning')
chk = DumbCheckpoint(file_dir+'prepare_continuous/manning', mode=FILE_READ)
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
    x_0, y_0, utm_zone, zone_letter = prepare_continuous.utm.from_latlon(lat, 0)
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
        200: {'elev': tidal_elev,'uv':tidal_v},  #set open boundaries to tidal_elev function
  }

def update_forcings(t):
  with timed_stage('update forcings'):
    print_output("Updating tidal field at t={}".format(t))
    elev = prepare_continuous.myboundary.set_tidal_field(Function(bathymetry2d.function_space()), t, dt)
    tidal_elev.project(elev) 
    v = prepare_continuous.myboundary.set_velocity_field(Function(VectorFunctionSpace(mesh2d,"CG",1)),t,dt)
    tidal_v.project(v)
    print_output("Done updating tidal field")

# a density function (to be optimised below) that specifies the number of turbines per unit area
turbine_density = Function(get_functionspace(mesh2d, "CG", 1), name='turbine_density')
# associate subdomain_id 2 (as in dx(2)) with a tidal turbine farm
# (implemented via a drag term) with specified turbine density
# Turbine characteristic can be specified via:
# - farm_options.turbine_options.thrust_coefficient (default 0.8)
# - farm_options.turbine_options.diameter (default 16.0)
farm_options = TidalTurbineFarmOptions()
farm_options.turbine_density = turbine_density
# amount of power produced per turbine (kW) on average to "break even" (cost = revenue)
# this is used to scale the cost, which is assumed to be linear with the number of turbines,
# in such a way that the cost is expressed in kW which can be subtracted from the profit
# which is calculated as the power extracted by the turbines
farm_options.break_even_wattage = BE/10
options.tidal_turbine_farms[turbine_area_PhyID] = farm_options

# we first run the "forward" model with no turbines
turbine_density.assign(0.0)



# create a density restricted to the farm
# the turbine_density, which is the control that will be varied in the optimisation,
# is defined everywhere, but it's influence is restricted to the farm area -
# the turbine drag term is integrated over the farm area only
# Because the turbine density is CG, the nodal values at the farm boundaries
# introduces a jagged edge around the farm where these values are tapered to zero.
# These nonzero values outside the farm itself do not contribute in any of the
# computations. For visualisation purposes we therefore project to a DG field
# restricted to the farm.
farm_density = Function(get_functionspace(mesh2d, "DG", 1), name='farm_density')
projector = SubdomainProjector(turbine_density, farm_density, turbine_area_PhyID)
projector.project()

cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')


# run as normal (this run will be annotated by firedrake_adjoint)
# solver_obj.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=Constant(0.0))
###spring:676,middle:492,neap:340###
solver_obj.load_state(492,outputdir='../../../outputs/0.validation/continuous-4cores')

solver_obj.iterate(update_forcings=update_forcings)


# compute maximum turbine density (based on a minimum of 2D between
# turbines laterally, and 5D in the streamwise direction)
D = farm_options.turbine_options.diameter
max_density = 1./(2*D*2*D)
print_output("Maximum turbine density = {}".format(max_density))

# we rescale the functional such that the gradients are ~ order magnitude 1.
# the scaling is based on the maximum cost term
# also we multiply by -1 so that if we minimize the functional, we maximize profit
# (maximize is also availble from pyadjoint but currently broken)
# scaling = -1/assemble(max(farm_options.break_even_wattage, 100) * max_density * dx(turbine_area_PhyID, domain=mesh2d))
# scaled_functional = scaling * cb.integrated_power
scaled_functional = -cb.average_profit[-1]
print_output(scaled_functional)
# specifies the control we want to vary in the optimisation
c = Control(turbine_density)

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
    optimisation.ControlsExportOptimisationCallback(solver_obj),
    optimisation.DerivativesExportOptimisationCallback(solver_obj),
    optimisation.UserExportOptimisationCallback(solver_obj, (farm_density,)),
    optimisation.FunctionalOptimisationCallback(solver_obj),
    turbines.TurbineOptimisationCallback(solver_obj, cb),
])

# this reduces the functional J(u, td) to a function purely of the control td:
# rf(td) = J(u(td), td) where the velocities u(td) of the entire simulation
# are computed by replaying the forward model for any provided turbine density td
rf = ReducedFunctional(scaled_functional, c, derivative_cb_post=callback_list)


if test_gradient:
    # whenever the forward model is changed - for example different terms in the equation,
    # different types of boundary conditions, etc. - it is a good idea to test whether the
    # gradient computed by the adjoint is still correct, as some steps in the model may
    # not have been annotated correctly. This can be done via the Taylor test.
    # Using the standard Taylor series, we should have (for a sufficiently smooth problem):
    #   rf(td0+h*dtd) - rf(td0) - < drf/dtd(rf0), h dtd> = O(h^2)

    direction_h= Function(turbine_density)
    direction_h.vector()[:]=numpy.random.random()*max_density

    # this tests whether the above Taylor series residual indeed converges to zero at 2nd order in h as h->0
    minconv = taylor_test(rf, turbine_density, direction_h)
    print_output("Order of convergence with taylor test (should be 2) = {}".format(minconv))

    assert minconv > 1.95

if optimise:
    # Optimise the control for minimal functional (i.e. maximum profit)
    # with a gradient based optimisation algorithm using the reduced functional
    # to replay the model, and computing its derivative via the adjoint
    # By default scipy's implementation of L-BFGS-B is used, see
    #   https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html
    # options, such as maxiter and pgtol can be passed on.
    chk = DumbCheckpoint(output_dir+'/optimal_density', mode=FILE_CREATE)
    td_opt = minimize(rf, bounds=[0, max_density],
                      options={'maxiter': 500, 'gtol': 1e-4})

    chk.store(td_opt, name='optimal_density')
    File(output_dir+'/optimal_density.pvd').write(td_opt)

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

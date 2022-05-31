"""
Using large domain's elevation and velocity outputs as boundary conditions.
"""

from thetis import *
import myboundary
import numpy
op2.init(log_level=INFO)

# setup the Thetis solver obj as usual:
mesh2d = Mesh('small_domain.msh')

tidal_amplitude = 5.
tidal_period = 12.42*60*60
H = 40


dt = 100.
t_export = 100
t_end = 20*100




# create solver and set options
solver_obj = solver2d.FlowSolver2d(mesh2d, Constant(H))
options = solver_obj.options
options.timestep = dt
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.output_directory = 'outputs'
options.check_volume_conservation_2d = True
options.element_family = 'dg-dg'
options.timestepper_type = 'CrankNicolson'
options.timestepper_options.implicitness_theta = 0.6
# using direct solver as PressurePicard does not work with dolfin-adjoint (due to .split() not being annotated correctly)
options.timestepper_options.solver_parameters = {'snes_monitor': None,
                                                 'snes_rtol': 1e-9,
                                                 'ksp_type': 'preonly',
                                                 'pc_type': 'lu',
                                                 'pc_factor_mat_solver_type': 'mumps',
                                                 'mat_type': 'aij'
                                                 }
options.horizontal_viscosity = Constant(100.0)
options.quadratic_drag_coefficient = Constant(0.0025)

# assign boundary conditions
left_tag = 1
right_tag = 2
coasts_tag = 3
tidal_elev = Function(FunctionSpace(mesh2d,"CG",1))
tidal_v = Function(VectorFunctionSpace(mesh2d,"CG",1))

solver_obj.bnd_functions['shallow_water'] = {
    left_tag: {'elev': tidal_elev, 'uv':tidal_v},
    right_tag: {'elev': tidal_elev, 'uv':tidal_v},
    coasts_tag: {'elev': tidal_elev, 'uv':tidal_v}
}

# a function to update the tidal_elev bc value every timestep
x = SpatialCoordinate(mesh2d)
g = 9.81
omega = 2 * pi / tidal_period


def update_forcings(t):
  with timed_stage('update forcings'):
    print_output("Updating tidal field at t={}".format(t))
    elev = myboundary.set_tidal_field(Function(FunctionSpace(mesh2d,"CG",1)), t, dt)
    tidal_elev.project(elev) 
    v = myboundary.set_velocity_field(Function(VectorFunctionSpace(mesh2d,"CG",1)), t, dt)
    tidal_v.project(v)
    print_output("Done updating tidal field")

solver_obj.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=tidal_elev)
locations,names = [(10000,2000)], ['point']
cb = DetectorsCallback(solver_obj, locations, ['elev_2d', 'uv_2d'], name='detectors',detector_names=names)
solver_obj.add_callback(cb, 'timestep')
solver_obj.iterate(update_forcings=update_forcings)


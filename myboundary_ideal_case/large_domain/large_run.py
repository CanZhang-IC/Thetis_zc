"""
Using large domain's elevation and velocity outputs as boundary conditions.
"""

from thetis import *
import numpy
op2.init(log_level=INFO)

# setup the Thetis solver obj as usual:
mesh2d = Mesh('large_domain.msh')

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
options.fields_to_export = ['uv_2d', 'elev_2d']
options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']
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
tidal_elev = Function(get_functionspace(mesh2d, "CG", 1), name='tidal_elev')
tidal_elev_bc = {'elev': tidal_elev}
# noslip currently doesn't work (vector Constants are broken in firedrake_adjoint)
freeslip_bc = {'un': Constant(0.0)}
solver_obj.bnd_functions['shallow_water'] = {
    left_tag: tidal_elev_bc,
    right_tag: tidal_elev_bc,
    coasts_tag: freeslip_bc
}

# a function to update the tidal_elev bc value every timestep
x = SpatialCoordinate(mesh2d)
g = 9.81
omega = 2 * pi / tidal_period


def update_forcings(t):
    print_output("Updating tidal elevation at t = {}".format(t))
    tidal_elev.project(tidal_amplitude*sin(omega*t + omega/pow(g*H, 0.5)*x[0]))

solver_obj.assign_initial_conditions(uv=as_vector((1e-7, 0.0)), elev=tidal_elev)
locations,names = [(10000,2000)], ['point']
cb = DetectorsCallback(solver_obj, locations, ['elev_2d', 'uv_2d'], name='detectors',detector_names=names)
solver_obj.add_callback(cb, 'timestep')
solver_obj.iterate(update_forcings=update_forcings)


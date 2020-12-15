"""
Idealised channel flow in 3D
============================

Solves shallow water equations in rectangular domain
with sloping bathymetry.

Flow is forced with tidal volume flux in the deep (ocean) end of the
channel, and a constant volume flux in the shallow (river) end.

This example demonstrates how to set up time dependent boundary conditions.
"""
from thetis import *
import time

t_start = time.time()

n_layers = 4
outputdir = '../../outputs/3D/dx150'
lx = 500
ly = 100
nx = 150
ny = 50
mesh2d = RectangleMesh(nx, ny, lx, ly)
print_output('Exporting to ' + outputdir)
t_end = 3600
t_export = 20

H = Constant(40)  # depth
P1_2d = FunctionSpace(mesh2d, "CG", 1)
bathymetry_2d = Function(P1_2d)
bathymetry_2d.assign(H)

u_in = 3.0
w_max = 5e-3

# create solver
solver_obj = solver.FlowSolver(mesh2d, bathymetry_2d, n_layers)
options = solver_obj.options
options.element_family = 'dg-dg'
options.timestepper_type = 'SSPRK22'
options.solve_salinity = False
options.solve_temperature = False
options.use_implicit_vertical_diffusion = True
options.use_bottom_friction = True
options.quadratic_drag_coefficient = Constant(0.0025)
options.use_ale_moving_mesh = False
options.use_lax_friedrichs_velocity = False
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.output_directory = outputdir
options.horizontal_velocity_scale = Constant(u_in)
options.vertical_velocity_scale = Constant(w_max)
options.fields_to_export = ['uv_2d', 'elev_2d', 'uv_3d',
                            'w_3d', 'uv_dav_2d']

solver_obj.create_fields()

xyz = SpatialCoordinate(solver_obj.mesh)
v_b = 50
v_inner = 0.1
v_length = 5
h_viscosity = Function(solver_obj.function_spaces.P1, name='viscosity')
h_viscosity.interpolate(conditional(ge(xyz[0],lx),(xyz[0]-(lx-v_length))*v_b/v_length+v_inner,v_inner))
File(outputdir+'/viscosity.pvd').write(h_viscosity)
options.horizontal_viscosity = h_viscosity
options.vertical_viscosity = h_viscosity

turbine_xyz = [50, 50, -20]
D = 20
th = 6  # thickness of disc
turbine_dims = [th, D, D]
C_T = 0.8
A_T = pi*(D/2)**2

bump = lambda x: exp(-1/Max(1-x**2, 0))
# bumps = [bump((x-turbine_x)/(dim/2)) for x, turbine_x, dim in zip(xyz, turbine_xyz, turbine_dims)]
# kernel = bumps[0]*bumps[1]*bumps[2]
radius_yz = sqrt((xyz[1]-turbine_xyz[1])**2 + (xyz[2]-turbine_xyz[2])**2)
kernel = bump((xyz[0] - turbine_xyz[0])/th)*bump(radius_yz/(D/2))

norm = assemble(kernel*dx)
kernel = kernel/Constant(norm)


u_inf = (solver_obj.fields.uv_3d + solver_obj.fields.uv_dav_3d) / Constant(0.5 * (1+sqrt(1-C_T)))
thrust_force = -Constant(0.5*C_T*A_T) * sqrt(dot(u_inf, u_inf)) * u_inf * kernel
# n = as_vector((sin(pi/6),cos(pi/6),0))
# thrust_force = -Constant(0.5*C_T*A_T) * abs(dot(u_inf, n)) * dot(u_inf,n) *n* kernel
options.momentum_source_3d = thrust_force

# boundary conditions are defined with a dict
# bnd conditions are assigned to each boundary tag with another dict
inflow_tag = 1
outflow_tag = 2
# assigning conditions for each equation
# these must be assigned before equations are created
u_ramped = Constant(0.0)
def update_forcings(t):
    u_ramped.assign(tanh(t/100.)*u_in)

solver_obj.bnd_functions['shallow_water'] = {inflow_tag: {'un': -u_ramped}, outflow_tag: {'elev': 0, 'un': u_ramped}}
solver_obj.bnd_functions['momentum'] = {inflow_tag: {'symm': None}, outflow_tag: {'symm': None}}

locations, names = [(70,50,-20)], ['1D']
cb = DetectorsCallback(solver_obj, locations, ['uv_3d','uv_dav_3d','w_3d'], name='detectors',detector_names=names)
solver_obj.add_callback(cb,'timestep')

solver_obj.iterate(update_forcings=update_forcings)

t_end = time.time()
print('time cost: {0:.2f}min'.format((t_end - t_start)/60))

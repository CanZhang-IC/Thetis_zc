from thetis import *

outputdir = 'nov11/outputs_v0.01'
mesh2d = Mesh('rectangular.msh')
print_output('Loaded mesh ' + mesh2d.name)
print_output('Exporting to ' + outputdir)

t_end = 500*20
t_export = 5*20

# bathymetry
P1_2d = FunctionSpace(mesh2d, 'CG', 1)
bathymetry_2d = Function(P1_2d, name='Bathymetry')
turbine_density = Function(FunctionSpace(mesh2d, "CG", 1), name='turbine_density').assign(0.0)

# Add viscosity sponge (depending on condition)
x = SpatialCoordinate(mesh2d)
h_viscosity = Function(P1_2d).interpolate(conditional(le(x[0], 3), 3.01-x[0], 0.01))

bathymetry_2d=Constant(0.54)
# --- create solver ---
solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry_2d)
options = solver_obj.options
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.output_directory = outputdir
options.check_volume_conservation_2d = True
options.fields_to_export = ['uv_2d', 'elev_2d','bathymetry_2d']
options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']
options.quadratic_drag_coefficient = Constant(0.0025)
options.timestepper_type = 'CrankNicolson'
options.horizontal_viscosity = h_viscosity
options.use_wetting_and_drying = True
options.wetting_and_drying_alpha = Constant(0.5)
if not hasattr(options.timestepper_options, 'use_automatic_timestep'):
    options.timestep = 20

options.timestepper_options.solver_parameters = {'snes_monitor': None,
                                                 'snes_rtol': 1e-5,
                                                 'ksp_type': 'preonly',
                                                 'pc_type': 'lu',
                                                 'pc_factor_mat_solver_type': 'mumps',
                                                 'mat_type': 'aij'
                                                 }

# Boundary conditions - Steady state case
tidal_elev = Function(P1_2d).assign(0.0)
tidal_vel = Function(P1_2d).assign(0.0)
solver_obj.bnd_functions['shallow_water'] = {1: {'un': tidal_vel},
                                             3: {'elev':tidal_elev}}

# initial conditions, piecewise linear function
elev_init = Function(P1_2d)
elev_init.assign(0.0)

def update_forcings(t_new):
    ramp = tanh(t_new / 2000.)
    tidal_vel.project(Constant(-ramp * 0.335))

# Initialise Discrete turbine farm characteristics
farm_options = DiscreteTidalTurbineFarmOptions()
farm_options.turbine_type = 'constant'
farm_options.turbine_options.thrust_coefficient = 0.6
farm_options.turbine_options.diameter = 0.27
farm_options.upwind_correction = False

farm_options.turbine_coordinates =[[Constant(25),Constant(0.6)]]
#add turbines to SW_equations
options.discrete_tidal_turbine_farms[1] = farm_options

#set initial condition
solver_obj.assign_initial_conditions(elev=elev_init, uv=(as_vector((1e-7, 0.0))))

# Operation of tidal turbine farm through a callback
cb = turbines.TurbineFunctionalCallback(solver_obj)
solver_obj.add_callback(cb, 'timestep')
locations=[]
names=[]
xx=[1,1.5,2,2.5,3,4,5,6,7,8]
yy=np.linspace(-1,1,300)
for i in xx:
    for j in yy:
        locations.append((25+0.27*i,0.6+0.27*j))
        names.append('name'+str((i,j)))  
cb2 = DetectorsCallback(solver_obj, locations, ['elev_2d', 'uv_2d'], name='detectors',detector_names=names)
solver_obj.add_callback(cb2, 'timestep')

# No update_forcings for steady state case
solver_obj.iterate(update_forcings=update_forcings)

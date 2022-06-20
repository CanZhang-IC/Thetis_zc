# calculation of 'bathymetry','distance for viscosity','manning','viscosity'
from thetis import * 
import utm
from scipy.io import netcdf_file
import scipy.interpolate
import numpy


mesh2d = Mesh('../mesh2/mesh.msh')


minimum_depth = 2 

utm_zone=51
utm_band='R'

def xy_to_latlon(xy, utm_zone,utm_band):
  '''
  convert the xy in utm into the xy in latlon
  '''
  latlon = []
  for x,y in xy:
    # when your mesh covers two utm_zones, make 'strict=False'
    latlon.append(utm.to_latlon(x, y, utm_zone, utm_band))
  return latlon

def get_bathymetry(bathymetry_file,bathymetry_file2, mesh2d):
  nc = netcdf_file(bathymetry_file)
  lat = nc.variables['lat'][:]
  lon = nc.variables['lon'][:]
  values = nc.variables['elevation'][:,:]
  #values = values.filled(9999.)
  interpolator = scipy.interpolate.RegularGridInterpolator((lat, lon), values)

  nc2 = netcdf_file(bathymetry_file2)
  lat2 = nc2.variables['lat'][:]
  lon2 = nc2.variables['lon'][:]
  values2 = nc2.variables['z'][:,:]
  #values = values.filled(9999.)
  interpolator2 = scipy.interpolate.RegularGridInterpolator((lat2, lon2), values2)

  P1_2d = FunctionSpace(mesh2d, 'CG', 1)
  bathymetry2d = Function(P1_2d, name="bathymetry")
  xvector = mesh2d.coordinates.dat.data
  bvector = bathymetry2d.dat.data
  assert xvector.shape[0]==bvector.shape[0]
  for i,xy in enumerate(xvector):
      lat, lon = utm.to_latlon(xy[0], xy[1], utm_zone, utm_band)
      
      
      if lat <= lat2.max() and lat >= lat2.min() and lon <= lon2.max() and lon >= lon2.min():
        bvector[i] = max(-interpolator2((lat, lon)), minimum_depth)
      else:
        bvector[i] = max(-interpolator((lat, lon)) ,minimum_depth)
  return bathymetry2d


def smoothen_bathymetry(bathymetry2d):  # smoothing bathymetry
  v = TestFunction(bathymetry2d.function_space())

  massb = assemble(v * bathymetry2d *dx)
  massl = assemble(v*dx)
  with massl.dat.vec as ml, massb.dat.vec as mb, bathymetry2d.dat.vec as sb:
      ml.reciprocal()
      sb.pointwiseMult(ml, mb)

chk = DumbCheckpoint('bathymetry', mode=FILE_CREATE)
with timed_stage('initialising bathymetry'):
  #bathymetry2d =get_bathymetry(main_dir+'/Netcdf/bathymetry.nc', mesh2d, 'z')
  bathymetry2d =get_bathymetry('../../Netcdf/bathymetry.nc', '../../Netcdf/bathymetry2.nc',mesh2d)

  smoothen_bathymetry(bathymetry2d)
  smoothen_bathymetry(bathymetry2d)
  chk.store(bathymetry2d, name='bathymetry')
  File('bathymetry.pvd').write(bathymetry2d)

###suitable bathymetry for turbine
def get_breakeven_bathymetry(max_depth,best_depth,min_depth,forbidden_depth):
  P1 = FunctionSpace(mesh2d,"CG",1)
  breakeven_bathymetry = Function(P1,name="breakeven_bathymetry")
  bvector = bathymetry2d.dat.data
  breakeven_bvector = breakeven_bathymetry.dat.data
  for i,depth in enumerate(bvector):
    if depth < forbidden_depth:
      breakeven_bvector[i] = 100
    elif depth < min_depth:
      breakeven_bvector[i] = 1+(min_depth-depth)/best_depth
    elif depth > max_depth:
      breakeven_bvector[i] = 1+(depth-max_depth)/best_depth
    else:
      breakeven_bvector[i] = 1
  return breakeven_bathymetry

chk = DumbCheckpoint('breakeven_bathymetry', mode=FILE_CREATE)
with timed_stage('initialising bathymetry'):
  max_depth,best_depth,min_depth,forbidden_depth = 30,30,30,20
  breakeven_bathymetry = get_breakeven_bathymetry(max_depth,best_depth,min_depth,forbidden_depth)
  smoothen_bathymetry(breakeven_bathymetry)
  smoothen_bathymetry(breakeven_bathymetry)
  chk.store(breakeven_bathymetry, name='breakeven_bathymetry')
  File('breakeven_bathymetry.pvd').write(breakeven_bathymetry)


# typical length scale
L = 1e3
V = FunctionSpace(mesh2d, 'CG', 1)

# Calculate distance to open boundary
print('Calculate distance for viscosity')

bcs = [DirichletBC(V, 0.0, 1000)] #make sure this matches physicalID of open boundaries

v = TestFunction(V)
u = Function(V)

solver_parameters = {
    'snes_type': 'ksponly',
    'ksp_rtol': 1e-4,
    'ksp_type': 'preonly',
    'pc_type': 'lu',
    'pc_factor_mat_solver_packages': 'mumps',
    # 'ksp_monitor_true_residual': True
    }

# Before we solve the Eikonal equation, let's solve a Laplace equation to
# generate an initial guess
F = L**2*(inner(grad(u), grad(v))) * dx - v * dx
solve(F == 0, u, bcs, solver_parameters=solver_parameters)

solver_parameters = {
     'snes_type': 'newtonls',
     'ksp_rtol': 1e-7,
         'ksp_type': 'preonly',
         'pc_type': 'lu',
         'pc_factor_mat_solver_packages': 'mumps',
         }

#epss values set the accuracy (in meters) of the final "distance to boundary" function. 
#To make more accurate add in extra iterations, eg, 500., 250., etc. This may result in the solver not converging.

epss = [100000., 10000., 5000., 2500., 1500., 1000.] 
# solve Eikonal equations
for i, eps in enumerate(epss):
  print('Solving Eikonal with eps == ' + str(float(eps)))
  F = inner(sqrt(inner(grad(u), grad(u))), v) * dx - v * dx + eps*inner(grad(u), grad(v)) * dx
  solve(F == 0, u, bcs, solver_parameters=solver_parameters)


chk = DumbCheckpoint('dist', mode=FILE_CREATE)
with timed_stage('initialising dist'):
  dist = Function(V, name='dist')
  dist.interpolate(u)
  chk.store(dist, name='dist')
  File('dist.pvd').write(u)

chk = DumbCheckpoint('manning', mode=FILE_CREATE)
with timed_stage('initialising manning'):
  manning = Function(V, name='manning')
  manning.interpolate(Max(0.02,0.1 * (1. - u / 1e2)))
  chk.store(manning, name='manning')
  File('manning.pvd').write(manning)

chk = DumbCheckpoint('viscosity', mode=FILE_CREATE)
with timed_stage('initialising viscosity'):
  h_viscosity = Function(V, name='viscosity')
  h_viscosity.interpolate(Max(1, 1000 * (1. - u / 1e2)))
  chk.store(h_viscosity, name='viscosity')
  File('viscosity.pvd').write(h_viscosity)

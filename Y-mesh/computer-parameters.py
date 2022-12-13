# calculation of 'bathymetry','distance for viscosity','manning','viscosity'
from cmath import nan
from thetis import * 
import utm
from scipy.io.netcdf import NetCDFFile
import scipy.interpolate
import numpy


mesh2d = Mesh('../Y-mesh/mesh.msh')


minimum_depth = -2 # 为正数时，表示不启用干湿边界条件。

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
  nc = NetCDFFile(bathymetry_file)
  lat = nc.variables['lat'][:]
  lon = nc.variables['lon'][:]
  values = nc.variables['elevation'][:,:]
  #values = values.filled(9999.)
  interpolator = scipy.interpolate.RegularGridInterpolator((lat, lon), values)

  nc2 = NetCDFFile(bathymetry_file2)
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

  massb = assemble(v * bathymetry2d *dx)#装配起来
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
      breakeven_bvector[i] = nan
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



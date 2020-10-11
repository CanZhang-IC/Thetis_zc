import uptide
import datetime
import utm
from math import tanh,sin,pi

#constituents = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', 'M4', 'MS4', 'MN4' ]
constituents = ['M2', 'S2']#'Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2']
tide = uptide.Tides(constituents)
tide.set_initial_time(datetime.datetime(2013,8,13,0,0,0)) #year, month, day, hour, min, sec

grid_file_name = "../Netcdf/grid_file.nc"
data_file_name = "../Netcdf/data_file.nc"
ranges=((120, 123),(29, 31))
tnci = uptide.tidal_netcdf.OTPSncTidalInterpolator(tide, grid_file_name, data_file_name, ranges)

utm_zone=51
utm_band='R'

def set_tidal_field(elev, t):
  tnci.set_time(t)
  mesh2d = elev.function_space().mesh()
  xvector = mesh2d.coordinates.dat.data# we have the coordinates of mesh
  evector = elev.dat.data# we want to use the coordinates to get the values on corresponding points 
  for i,xy in enumerate(xvector):
    lat, lon = utm.to_latlon(xy[0], xy[1], utm_zone, utm_band)
    try:
      evector[i] = tnci.get_val((lon, lat))
    except uptide.netcdf_reader.CoordinateError:
      evector[i] = 0.
  return elev

def set_tidal_field1(elev, t):
  tnci.set_time(t)
  mesh2d = elev.function_space().mesh()
  xvector = mesh2d.coordinates.dat.data# we have the coordinates of mesh
  evector = elev.dat.data# we want to use the coordinates to get the values on corresponding points 
  for i,xy in enumerate(xvector):
    evector[i] = sin(2 * pi /12.42*60*60*t + 2 * pi /12.42*60*60/pow(9.81*40, 0.5)*xy[0])

def set_tidal_field_ramp(elev, t, start):
    tnci.set_time(t)
    mesh2d = elev.function_space().mesh()
    xvector = mesh2d.coordinates.dat.data
    evector = elev.dat.data
    ramp = tanh((t-start)/5000.)
    for i, xy in enumerate(xvector):
        lat, lon = utm.to_latlon(xy[0], xy[1], utm_zone, utm_band)
        try:
            evector[i] = tnci.get_val((lon, lat))*ramp    # Adding initial a correction depth for LAT
        except uptide.netcdf_reader.CoordinateError:
            evector[i] = 0.    # Adding initial a correction depth for LAT

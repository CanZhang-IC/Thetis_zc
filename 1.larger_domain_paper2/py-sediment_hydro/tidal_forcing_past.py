import uptide
import uptide.tidal_netcdf
import datetime
import utm

main_dir = '../../..'
#constituents = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', 'M4', 'MS4', 'MN4' ]
constituents = ['Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2']
tide = uptide.Tides(constituents)
tide.set_initial_time(datetime.datetime(2013,8,8,0,0,0)) #year, month, day, hour, min, sec

grid_file_name = main_dir+"/Netcdf/grid_file.nc"
data_file_name = main_dir+"/Netcdf/data_file.nc"
ranges=((121, 124),(28, 31)) #lon/lat 
tnci = uptide.tidal_netcdf.OTPSncTidalInterpolator(tide, grid_file_name, data_file_name, ranges)

utm_zone=51
utm_band='R'

def set_tidal_field(elev, t):
  tnci.set_time(t)
  mesh2d = elev.function_space().mesh()
  xvector = mesh2d.coordinates.dat.data# we have the coordinates of mesh
  evector = elev.dat.data# we want to use the coordinates to get the values on corresponding points 
  for i,xy in enumerate(xvector):

    # outproj = pyproj.Proj(init='epsg:4326') # needs to match your qgis layers (bottom right hand corner of the qgis display window)
    # outproj = pyproj.Proj(init='epsg:32651') # needs to be aligned with your utm zone
    # inproj = pyproj.Proj(init='epsg:32651') # needs to be aligned with your utm zone
    # lon, lat =pyproj.transform(inproj, outproj, xy[0], xy[1])
    lat, lon = utm.to_latlon(xy[0], xy[1], utm_zone, utm_band)
    #print('lat:',lat,'lon:',lon)
    try:
      evector[i] = tnci.get_val((lon, lat))
    except uptide.netcdf_reader.CoordinateError:
      evector[i] = 0.
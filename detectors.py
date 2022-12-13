from thetis import *
import pyproj
import numpy as np

tidegauge_file = 'tidal_gauges.csv'  # name of your file

UTM_ZONE = pyproj.Proj(
        proj='utm',
        zone=51,  # utm zone number
        datum='WGS84',
        units='m',
        errcheck=True)
LL_WGS84 = pyproj.Proj(proj='latlong', datum='WGS84', errcheck=True)

def get_detectors(mesh2d):
    gauge_names = np.loadtxt(tidegauge_file, skiprows=1, usecols=(0,), dtype=str, delimiter=',')
    gauge_latlon = np.loadtxt(tidegauge_file, skiprows=1, usecols=(1,2), delimiter=',') #your lat, lon data is in cols 2/3
    ind = np.argsort(gauge_names)#Returns the indices that would sort an array.
    gauge_names = list(gauge_names[ind])
    gauge_latlon = list(gauge_latlon[ind])

    #make names unique
    unique_names = []
    last_name = ''; ctr = 0
    for name in gauge_names:
        if name==last_name:
            unique_names.append(name + '_' + str(ctr))
            ctr += 1
        else:
            unique_names.append(name)
            ctr = 0
        last_name = name

    inproj = pyproj.Proj(init='epsg:4326') #lat / lon coords
    outproj = pyproj.Proj(init='epsg:32651') # 326 + utm zone number (51 for Korea)
    gauge_xy = [pyproj.transform(inproj, outproj, lon, lat) for lat, lon in gauge_latlon] # convert latlon to utm xy
    
    return select_and_move_detectors(mesh2d, gauge_xy, unique_names, maximum_distance=10e3)#Select those detectors that are within the domain and/or move them to the nearest cell centre within the domain

if __name__ == "__main__":
    inproj = pyproj.Proj(init='epsg:4326') #lat / lon coords
    outproj = pyproj.Proj(init='epsg:32651') # 326 + utm zone number (51 for Korea)
    gauge_latlon = [(30.0334,122.4098),(30.0373,122.4149)]
    gauge_xy = [pyproj.transform(inproj, outproj, lon, lat) for lat, lon in gauge_latlon]
    print(gauge_xy)
    # mesh2d = Mesh('./mesh/mesh.msh')  # mesh file
    # # import pandas as pd
    # # df_file = pd.read_excel('./location.xlsx',usecols="M:N")
    # # print(df_file.keys())
    # # xx = np.array(df_file['X.1'][89:107])
    # # yy = np.array(df_file['Y.1'][89:107])
    # xmin,ymin,xmax,ymax = 443340, 3322634, 443592, 3322848 
    # locations = [[444000,3323000]]
    # for x in range(xmin+20,xmax-20,60):
    #     for y in range(ymin+20,ymax-20,120):
    #         locations.append([x,y])
    # for x in range(xmin+20+30,xmax-20,60):
    #     for y in range(ymin+20+60,ymax-20,120):
    #         locations.append([x,y])
    # names = []
    # for i in range(len(locations)):
    #     names.append('point ' + str(i))
    # # locations = []
    # # for x,y in zip(xx,yy):
    # #     locations.append([x,y])
    # # print(list(locations))


    # # locations, names = get_detectors(mesh2d)
    # if mesh2d.comm.rank == 0: # only processor 0
    #     print_output("Found detectors: {}".format(names))
    #     # write out shape-file
    #     import shapely.geometry
    #     import fiona
    #     import fiona.crs

    #     schema = {'geometry': 'Point', 'properties': {'name': 'str'}}
    #     crs = fiona.crs.from_string(UTM_ZONE.srs)
    #     shpfile= fiona.collection("detectors.shp", "w", "ESRI Shapefile", schema, crs=crs) 
    #     for xy, name in zip(locations, names):
    #         print(xy,name)
    #         point = shapely.geometry.Point(xy[0], xy[1])
    #         shpfile.write({'properties': {'name': 'point'}, 'geometry': shapely.geometry.mapping(point)})


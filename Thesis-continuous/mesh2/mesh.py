import numpy as np
import qmesh

qmesh.initialise()
#Reading in the shapefile describing the domain boundaries, and creating a gmsh file.
boundaries = qmesh.vector.Shapes()
boundaries.fromFile('turbcircle.shp')  # Domain file
loopShapes = qmesh.vector.identifyLoops(boundaries,isGlobal=False, defaultPhysID=1000,fixOpenLoops=True)
polygonShapes = qmesh.vector.identifyPolygons(loopShapes, smallestNotMeshedArea=10, meshedAreaPhysID = 1)

innerSound_plot_lines = qmesh.vector.Shapes()
innerSound_plot_lines.fromFile('tl.shp')
innerSound_plot_loops = qmesh.vector.identifyLoops(innerSound_plot_lines,isGlobal=False, defaultPhysID=500,
          fixOpenLoops=True, extraPointsPerVertex=10)
innerSound_plot_polygon = qmesh.vector.identifyPolygons(innerSound_plot_loops, smallestNotMeshedArea=50,meshedAreaPhysID = 2)

#Insert tidal plots and turbines into domain.
domainLines, domainPolygons = qmesh.vector.insertRegions(
                                      loopShapes, polygonShapes,\
                                      innerSound_plot_loops, innerSound_plot_polygon)

#Create raster for mesh gradation 
Turbine_boundaries = qmesh.vector.Shapes()
Turbine_boundaries.fromFile('turbturb.shp') 
gradationRaster_shoreline_1 = qmesh.raster.gradationToShapes()
gradationRaster_shoreline_1.setShapes(Turbine_boundaries)
gradationRaster_shoreline_1.setRasterBounds(122.0,122.8,29.7,30.3)  # lonmin, lonmax, latmin, latmax
gradationRaster_shoreline_1.setRasterResolution(1500,1500)
gradationRaster_shoreline_1.setGradationParameters(100.0,1500.0,0.1,0.001)  # resmin, resmax (m), rate of change
gradationRaster_shoreline_1.calculateLinearGradation()
gradationRaster_shoreline_1.writeNetCDF('gradation_to_turbines_lines.nc')

#Create raster for mesh gradation towards Inner Sound
gradationRaster_shoreline_2 = qmesh.raster.gradationToShapes()
gradationRaster_shoreline_2.setShapes(innerSound_plot_polygon)
gradationRaster_shoreline_2.setRasterBounds(122.0,122.8,29.7,30.3)
gradationRaster_shoreline_2.setRasterResolution(1500,1500)
gradationRaster_shoreline_2.setGradationParameters(100,1500.0,0.1) 
gradationRaster_shoreline_2.calculateLinearGradation()


#meshMetricRaster.writeNetCDF('meshMetric.nc)
# meshMetricRaster = gradationRaster_shoreline_1
meshMetricRaster = qmesh.raster.meshMetricTools.minimumRaster([gradationRaster_shoreline_1, gradationRaster_shoreline_2])
domain = qmesh.mesh.Domain()
# domain.setGeometry(loopShapes, polygonShapes)
domain.setGeometry(domainLines, domainPolygons)
domain.setMeshMetricField(meshMetricRaster)
domain.setTargetCoordRefSystem('EPSG:32651', fldFillValue=1000.0)
domain.gmsh(geoFilename='mesh.geo', \
            fldFilename='mesh.fld', \
            mshFilename='mesh.msh')

mesh = qmesh.mesh.Mesh()
mesh.readGmsh('mesh.msh', 'EPSG:32651')  # 326 + UTM number
mesh.writeShapefile('mesh')


import numpy as np
import qmesh

qmesh.initialise()
#Reading in the shapefile describing the domain boundaries, and creating a gmsh file.
boundaries = qmesh.vector.Shapes()
boundaries.fromFile('largelines.shp')  # Domain file
loopShapes = qmesh.vector.identifyLoops(boundaries,isGlobal=False, defaultPhysID=1000,fixOpenLoops=True)
polygonShapes = qmesh.vector.identifyPolygons(loopShapes, smallestNotMeshedArea=10, 
                                                                 meshedAreaPhysID = 1)
polygonShapes.writeFile('polygonshape.shp')

domainLines, domainPolygons = loopShapes, polygonShapes

innerSound_plot_lines = qmesh.vector.Shapes()
innerSound_plot_lines.fromFile('tl1.shp')
innerSound_plot_loops = qmesh.vector.identifyLoops(innerSound_plot_lines,isGlobal=False, defaultPhysID=1000,
          fixOpenLoops=True, extraPointsPerVertex=10)
innerSound_plot_polygon = qmesh.vector.identifyPolygons(innerSound_plot_loops, smallestNotMeshedArea=50,
                                                                 meshedAreaPhysID = 2)
                                 

#Insert tidal plots and turbines into domain.
domainLines1, domainPolygons1 = qmesh.vector.insertRegions(
                                      loopShapes, polygonShapes,\
                                      innerSound_plot_loops, innerSound_plot_polygon)

innerSound_plot_lines2 = qmesh.vector.Shapes()
innerSound_plot_lines2.fromFile('tl2.shp')
innerSound_plot_loops2 = qmesh.vector.identifyLoops(innerSound_plot_lines2,isGlobal=False, defaultPhysID=1000,
          fixOpenLoops=True, extraPointsPerVertex=10)
innerSound_plot_polygon2 = qmesh.vector.identifyPolygons(innerSound_plot_loops2, smallestNotMeshedArea=50,
                                                                 meshedAreaPhysID = 2)


domainLines2, domainPolygons2 = qmesh.vector.insertRegions(
                                      domainLines1, domainPolygons1,\
                                      innerSound_plot_loops2, innerSound_plot_polygon2)

innerSound_plot_lines3 = qmesh.vector.Shapes()
innerSound_plot_lines3.fromFile('tl3.shp')
innerSound_plot_loops3 = qmesh.vector.identifyLoops(innerSound_plot_lines3,isGlobal=False, defaultPhysID=1000,
          fixOpenLoops=True, extraPointsPerVertex=10)
innerSound_plot_polygon3 = qmesh.vector.identifyPolygons(innerSound_plot_loops3, smallestNotMeshedArea=50,
                                                                 meshedAreaPhysID = 2)


domainLines3, domainPolygons3 = qmesh.vector.insertRegions(
                                      domainLines2, domainPolygons2,\
                                      innerSound_plot_loops3, innerSound_plot_polygon3)

innerSound_plot_lines4 = qmesh.vector.Shapes()
innerSound_plot_lines4.fromFile('tl4.shp')
innerSound_plot_loops4 = qmesh.vector.identifyLoops(innerSound_plot_lines4,isGlobal=False, defaultPhysID=1000,
          fixOpenLoops=True, extraPointsPerVertex=10)
innerSound_plot_polygon4 = qmesh.vector.identifyPolygons(innerSound_plot_loops4, smallestNotMeshedArea=50,
                                                                 meshedAreaPhysID = 3)


domainLines, domainPolygons = qmesh.vector.insertRegions(
                                      domainLines3, domainPolygons3,\
                                      innerSound_plot_loops4, innerSound_plot_polygon4)
#Create raster for mesh gradation 
GSHHS_fine_boundaries1 = qmesh.vector.Shapes()
GSHHS_fine_boundaries1.fromFile('largecoastlines1.shp')  # Coastline file
gradationRaster_shoreline1 = qmesh.raster.gradationToShapes()
gradationRaster_shoreline1.setShapes(GSHHS_fine_boundaries1)
gradationRaster_shoreline1.setRasterBounds(121.0, 123.4, 29.4, 30.9)  # lonmin, lonmax, latmin, latmax
gradationRaster_shoreline1.setRasterResolution(1500,1500)
gradationRaster_shoreline1.setGradationParameters(500.0,3000.0,0.5)  #184792 nodes, 195329 vertices 382186 elements
gradationRaster_shoreline1.calculateLinearGradation()


GSHHS_fine_boundaries = qmesh.vector.Shapes()
GSHHS_fine_boundaries.fromFile('largecoastlines2.shp')  # Coastline file
gradationRaster_shoreline = qmesh.raster.gradationToShapes()
gradationRaster_shoreline.setShapes(GSHHS_fine_boundaries)
gradationRaster_shoreline.setRasterBounds(121.0, 123.4, 29.4, 30.9)  # lonmin, lonmax, latmin, latmax
gradationRaster_shoreline.setRasterResolution(1500,1500)
gradationRaster_shoreline.setGradationParameters(300.0,3000.0,0.5)  #50213 nodes, 60750 vertices 112140 elements
gradationRaster_shoreline.calculateLinearGradation()


re = 200
#Create raster for mesh gradation towards Inner Sound
gradationRaster_shoreline_1 = qmesh.raster.gradationToShapes()
gradationRaster_shoreline_1.setShapes(innerSound_plot_polygon)
gradationRaster_shoreline_1.setRasterBounds(121.0, 123.4, 29.4, 30.9)
gradationRaster_shoreline_1.setRasterResolution(1500,1500)
gradationRaster_shoreline_1.setGradationParameters(re,3000.0,0.01) 
gradationRaster_shoreline_1.calculateLinearGradation()



gradationRaster_shoreline_2 = qmesh.raster.gradationToShapes()
gradationRaster_shoreline_2.setShapes(innerSound_plot_polygon2)
gradationRaster_shoreline_2.setRasterBounds(121.0, 123.4, 29.4, 30.9)
gradationRaster_shoreline_2.setRasterResolution(1500,1500)
gradationRaster_shoreline_2.setGradationParameters(re,3000.0,0.01) 
gradationRaster_shoreline_2.calculateLinearGradation()

gradationRaster_shoreline_3 = qmesh.raster.gradationToShapes()
gradationRaster_shoreline_3.setShapes(innerSound_plot_polygon3)
gradationRaster_shoreline_3.setRasterBounds(121.0, 123.4, 29.4, 30.9)
gradationRaster_shoreline_3.setRasterResolution(1500,1500)
gradationRaster_shoreline_3.setGradationParameters(re,3000.0,0.01) 
gradationRaster_shoreline_3.calculateLinearGradation()

gradationRaster_shoreline_4 = qmesh.raster.gradationToShapes()
gradationRaster_shoreline_4.setShapes(innerSound_plot_polygon4)
gradationRaster_shoreline_4.setRasterBounds(121.0, 123.4, 29.4, 30.9)
gradationRaster_shoreline_4.setRasterResolution(1500,1500)
gradationRaster_shoreline_4.setGradationParameters(re,3000.0,0.01) 
gradationRaster_shoreline_4.calculateLinearGradation()




#meshMetricRaster.writeNetCDF('meshMetric.nc)
#meshMetricRaster = gradationRaster_shoreline
meshMetricRaster = qmesh.raster.meshMetricTools.minimumRaster([gradationRaster_shoreline, gradationRaster_shoreline1,
    gradationRaster_shoreline_1,
	gradationRaster_shoreline_2,
    gradationRaster_shoreline_3,
    gradationRaster_shoreline_4,
    ])
domain = qmesh.mesh.Domain()
#domain.setGeometry(loopShapes, polygonShapes)
meshname='mesh'
domain.setGeometry(domainLines, domainPolygons)
domain.setMeshMetricField(meshMetricRaster)
domain.setTargetCoordRefSystem('EPSG:32651', fldFillValue=1000.0)
domain.gmsh(geoFilename=meshname +'.geo', \
            fldFilename=meshname +'.fld', \
            mshFilename=meshname +'.msh')

mesh = qmesh.mesh.Mesh()
mesh.readGmsh(meshname +'.msh', 'EPSG:32651')  # 326 + UTM number
mesh.writeShapefile(meshname)


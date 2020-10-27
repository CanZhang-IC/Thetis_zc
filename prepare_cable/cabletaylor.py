from thetis import *
from firedrake_adjoint import *
from cable_overloaded import cablelength
import numpy
import Hybrid_Code


turbine_locations = [[x,y] for x in numpy.arange(site_x1+20, site_x2-20, 60) for y in numpy.arange(site_y1+20, site_y2-20, 50)]
landpointlocation = [[444000,3323000]]
#cableclass = Hybrid_Code.CableCostGA(turbine_locations, substation_location=landpointlocation,show_prog = False, show_result = False)
order_w = [[0, 1, 6, 10, 11, 7, 8, 4], [0, 5, 9, 13, 14, 15, 16, 12], [0, 2, 3]]#cableclass.compute_cable_cost_order()
order_con = [Constant(i) for j in order_w for i in j]
turbine_locations_con = [Constant(x) for xy in turbine_locations for x in xy]
landpointlocation_con = [Constant(x) for xy in landpointlocation for x in xy]

J = cablelength(turbine_locations_con,landpointlocation_con,order_con)

c = [Control(x) for x in turbine_locations_con] 

rf = ReducedFunctional(J , c)



print(J, rf(turbine_locations_con)) # should be the same

for i,xy in enumerate(turbine_locations_con):
    xy.assign(Constant(i*100))



print(rf(turbine_locations_con))

print(cablelength(turbine_locations_con,landpointlocation_con,order_con))

# h0=[]
# for i in range(len(turbine_locations)):
#     h0.append(Constant(20))
#     h0.append(Constant(10))

# minconv = taylor_test(rf, turbine_locations_con, h0)
# print_output("Order of convergence with taylor test (should be 2) = {}".format(minconv))

# assert minconv > 1.95
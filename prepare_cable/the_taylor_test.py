from thetis import *
from firedrake_adjoint import *
from centretoland_overloaded import two_length
import numpy


turbine_locations = [[x,y] for x in numpy.arange(443032+20, 443288-20, 60) for y in numpy.arange(3322891+20, 3323091-20, 40)]
landpointlocation = [[442500,3322750]]


turbine_locations_con = [Constant(x) for xy in turbine_locations for x in xy]
landpointlocation_con = [Constant(x) for xy in landpointlocation for x in xy]

J = two_length(turbine_locations_con,landpointlocation_con)

c = [Control(x) for x in turbine_locations_con] 

rf = ReducedFunctional(J , c)



# print(J, rf(turbine_locations_con)) # should be the same

# for i,xy in enumerate(turbine_locations_con):
#     xy.assign(Constant(i*100))

# print(rf(turbine_locations_con))

# print(two_length(turbine_locations_con,landpointlocation_con))

h0=[]
for i in range(len(turbine_locations)):
    h0.append(Constant(10))
    h0.append(Constant(10))

minconv = taylor_test(rf, turbine_locations_con, h0)
print_output("Order of convergence with taylor test (should be 2) = {}".format(minconv))

assert minconv > 1.95
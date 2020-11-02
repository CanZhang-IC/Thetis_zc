"""
Discrete turbines optimisation example
=======================================
Test 1: Forward example
Test 2: Forward example (change interpolate to project)
Test 3: Taylor test
"""
# to enable a gradient-based optimisation using the adjoint to compute
# gradients, we need to import from thetis_adjoint instead of thetis. This
# ensure all firedrake operations in the Thetis model are annotated
# automatically, in such a way that we can rerun the model with different input
# parameters, and also derive the adjoint-based gradient of a specified input
# (the functional) with respect to a specified input (the control)
from thetis import *
from firedrake_adjoint import *
from pyadjoint import minimize
import numpy
op2.init(log_level=INFO)


mesh2d = Mesh('headland2.msh')
P1_2d = FunctionSpace(mesh2d, 'CG', 1)
angle_function = Function(P1_2d)
x = SpatialCoordinate(mesh2d)
radius = 20

positions = [953.357502692696, 280.9002564543021, 930.0000000000358, 320.6500260306605, 965.0082634158863, 340.0, 998.0989933583933, 267.5617248267398, 985.0082634158797, 305.3589838486847, 1005.0082634158774, 339.9999999999895, 1037.377744836081, 260.0, 1024.2870148935542, 297.7972590219382, 1044.7315589859368, 335.30321502951915]
turbine_coordinates =[[Constant(positions[i]), Constant(positions[i+1])] for i in range(0,len(positions),2)]

alphas = [Constant(i) for i in 
[0.0, 0.0, 0.23915608608376898, 3.925031370943139, 3.626463857267974, 10.642066063788867, 7.385338800655825, 8.466238843590158, 15.28695703243358]]

for coord, alpha in zip(turbine_coordinates,alphas):
    rho = sqrt(x[0]**2 + x[1]**2)
    theta = atan_2(x[0],x[1]) + alpha/180*pi
    rho2 = sqrt(coord[0]**2 + coord[1]**2)
    theta2 = atan_2(coord[0],coord[1]) + alpha/180*pi
    dx0 = (rho * cos(theta) - rho2 * cos(theta2))/radius
    dx1 = (rho * sin(theta) - rho2 * sin(theta2))/radius
    bump = conditional(lt(abs(dx0), 1), conditional(lt(abs(dx1), 0.5), exp(1-1/(1-dx1**2))*exp(1-1/(1-dx0**2)), 0), 0)
    angle_function.interpolate(angle_function+bump/1.09665)
    # dx0 = (x[0] - coord[0])/radius
    # dx1 = (x[1] - coord[1])/radius
    # dis_xy = dx0*dx0+dx1*dx1
    # bump = conditional(lt(dis_xy,1),exp(1-1/(1-dis_xy)), 0)
    # angle_function.interpolate(angle_function+bump/1.2681)

File('without_headland2.pvd').write(angle_function)


    

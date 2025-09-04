from thetis import *
from firedrake_adjoint import *
import scipy.interpolate
import numpy as np

#choosing a rectangular area which covers the small domain and inside the large domain.
#This is because the scipy.interpolate is not good at working with unstructured grid data.
#So we set a rectangular area with high resolution by firedrake.at function and
#then the scipy.interpolate can work well with the data inside the rectangular.
min_x, min_y,max_x,max_y= 428500,3306000,460000,3335511
n = 100 #ensure adequate mesh points to refine our solution
gridxy1 = np.mgrid[min_x:max_x:100j, min_y:max_y:100j].T
gridxy = np.reshape(gridxy1, (n*n, 2))
gridx, gridy = gridxy.T

###Read in the large domain's mesh
large_mesh2d = Mesh('../../1.larger_domain_paper2/mesh/mesh.msh')# this is where you should make changes 

###Functionspace for elevation
P1 = FunctionSpace(large_mesh2d, "DG", 1)
###Functionspace for velocity
VP = VectorFunctionSpace(large_mesh2d, "DG", 1)

#location where the elevation and velocity outputs from the large domain stores
h5file_dir = '../../../outputs/paper2validation-4cores/hdf5' # this is where you should make changes 

# Be careful here!!!
# The 'dt' in the small domain's script must be the same with 't_export' in the large domain's script
# This is to make sure the large domain's outputs can provide elevation and velocity results for every timestep for the small domain 

def set_tidal_field(elev, t, dt):
    num = str(int(t/dt))
    while len(num) < 5:
        num = '0' + num

    chk = DumbCheckpoint(h5file_dir+'/Elevation2d_'+num, mode=FILE_READ)
    ele_large = Function(P1)
    chk.load(ele_large, name='elev_2d')
    chk.close()

    ele_l_xy = ele_large.at(gridxy,dont_raise=True)
    ele_l_xy = [np.array([0]) if x == None else x for x in ele_l_xy]
    ele_l_xy = np.array(ele_l_xy, dtype=object)
    ele_l_xy.astype(np.float64)

    interpolator = scipy.interpolate.LinearNDInterpolator(gridxy, ele_l_xy)

    mesh2d = elev.function_space().mesh()
    xvector = mesh2d.coordinates.dat.data
    evector = elev.dat.data

    for i,xy in enumerate(xvector):
        if interpolator(xy) == None:
            evector[i] = 0
        else:
            evector[i] = interpolator(xy)
    return elev

def set_velocity_field(uv_2d,t,dt):
    num = str(int(t/dt))
    while len(num) < 5:
        num = '0' + num

    chk = DumbCheckpoint(h5file_dir+'/Velocity2d_'+num, mode=FILE_READ)
    v_large = Function(VP)
    chk.load(v_large, name='uv_2d')
    chk.close()

    v_l_xy = v_large.at(gridxy,dont_raise=True)
    

    v_l_x = [x[0] if isinstance(x,np.ndarray) else np.array([0]) for x in v_l_xy]
    v_l_x = np.array(v_l_x, dtype=object)
    v_l_x.astype(np.float64)

    v_l_y = [x[1] if isinstance(x,np.ndarray) else np.array([0]) for x in v_l_xy]
    v_l_y = np.array(v_l_y, dtype=object)
    v_l_y.astype(np.float64)

    xinterpolator = scipy.interpolate.LinearNDInterpolator(gridxy, v_l_x)
    yinterpolator = scipy.interpolate.LinearNDInterpolator(gridxy, v_l_y)

    mesh2d = uv_2d.function_space().mesh()
    xvector = mesh2d.coordinates.dat.data
    evector = uv_2d.dat.data

    for i,xy in enumerate(xvector):
        vx = xinterpolator(xy)
        vy = yinterpolator(xy)
        if vx == None:
            vx = 0
        if vy == None:
            vy = 0
        evector[i] = [vx,vy]
    return uv_2d


    


if __name__ == "__main__":
    pass

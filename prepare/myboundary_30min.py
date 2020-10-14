from thetis import *
from firedrake_adjoint import *
import scipy.interpolate
import numpy as np

min_x, min_y,max_x,max_y= 420000,3300000,470000,3340000
n = 100 #ensure adequate mesh points to refine our solution
gridxy1 = np.mgrid[min_x:max_x:100j, min_y:max_y:100j].T
gridxy = np.reshape(gridxy1, (n*n, 2))
gridx, gridy = gridxy.T

large_mesh2d = Mesh('/media/can/can_disk/thetis_new/paper2/continuemethod/mesh/mesh.msh')
P1 = FunctionSpace(large_mesh2d, "DG", 1)
VP = VectorFunctionSpace(large_mesh2d, "DG", 1)

h5file_dir = '/media/can/can_disk/thetis_new/paper2/continuemethod/outputs/paper2validation/hdf5'

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

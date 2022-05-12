import h5py
import matplotlib.pyplot as plt
import os
import numpy as np
from thetis import *
from firedrake_adjoint import *

mesh2d = Mesh('/home/can/Git_thetis/sediment_hydro/mesh.msh')
D = FunctionSpace(mesh2d, "CG", 1)
V = VectorFunctionSpace(mesh2d,"DG",1)

# path = os.walk('./')
# special_marks=[]
# for root,dirs,files in path:
#     special_marks.append(dirs)
special_marks = [['sediment_size-v1e-9']]
for special_mark in special_marks[0]:
    result_dir = './'+special_mark+'/'
    print(result_dir)
    filenames = range(4,41,4)

    plt.figure(figsize=(16,9))
    for pa in filenames:
        path = os.walk(result_dir+str(float(pa))+'/hdf5/')
        a=[]
        for root,dirs,files in path:
            a.append(files)
        iteration_number_all = []
        for name in a[0] :
            # print(name)
            iteration_number_all.append(int(name[-6:-3]))
        last_iteration = max(iteration_number_all)
        print(pa,'\t',last_iteration)
        
        if last_iteration < 45:
            label = str(pa/100)+'mm-crashed'
        else:
            label = str(pa/100)+'mm'
        # label = 'viscosity' + str(pa)
        while len(str(last_iteration)) < 5:
            last_iteration = '0' + str(last_iteration)
        chk = DumbCheckpoint(result_dir+str(float(pa))+'/hdf5/bathymetry2d_'+str(last_iteration),mode=FILE_READ)
        bathymetry_2d = Function(D, name='bathymetry_2d')
        chk.load(bathymetry_2d)
        chk.close()
        File(result_dir+'depth/depth.pvd',mode='a').write(bathymetry_2d)


        chk_v = DumbCheckpoint(result_dir+str(float(pa))+'/hdf5/Velocity2d_'+str(last_iteration),mode=FILE_READ)
        velocity_final = Function(V, name='uv_2d')
        chk_v.load(velocity_final)
        chk_v.close()
        File(result_dir+'velocity/velocity_final.pvd',mode='a').write(velocity_final)
        # record final bathymetry for plotting
        xaxisthetis1 = []
        baththetis1 = []

        for i in np.linspace(100, 950, 1000):
            xaxisthetis1.append(i)
            baththetis1.append(-bathymetry_2d.at([i, 150]))

        plt.plot(xaxisthetis1, baththetis1, label=label)

    plt.legend()
    plt.savefig(result_dir[:-2]+'-fig.png',dpi=300)

    # plt.figure()
    # for pa in filenames:
    #     def_file = h5py.File(result_dir+str(float(pa))+'/diagnostic_turbine.hdf5','r+')
    #     a = []
    #     for name,data in def_file.items():
            
    #         if name == 'current_power':
    #             print(data[0],data[-1])
    #             for i in data:
    #                 a.append(i)
    #     b = range(len(a))
    #     plt.plot(b,a,'.-',label=str(pa/100)+'mm')
    # plt.legend()
    # plt.savefig(result_dir+'power.png',dpi=300)


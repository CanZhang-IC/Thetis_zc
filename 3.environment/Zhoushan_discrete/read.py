import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
import os


thedir=os.getcwd()  

filenames = ['controls','velocity_errors','velocity_errors2','each_farm_optimisation','eachturbine','farm_optimisation','functional','turbine','volume2d']

BE = 2.0
for i in [1,2]:
    output_dir =  '../../../outputs/3.environment/discrete/flood_ebb-forward/y-BE-'+str(BE)[:-2]
    df = h5py.File(output_dir+'/diagnostic_'+filenames[i]+'.hdf5','r+')
    # df = h5py.File(output_dir+'/diagnostic_sediment_errors.hdf5','r+')

    for name,data in df.items():
        if name == 'time':
            pass
        else:
            print(name)

            print(data[0],'\n',len(data),data[-1])
            if name == 'RMSE_current':
                a = int(len(data)/6)
                b = np.linspace(492,491+a,a)
                plt.plot(b,data[::6],'.-',label='2op-'+str(i))
                plt.legend()

BE = 4.0
for i in [1,2]:
    output_dir =  '../../../outputs/3.environment/discrete/flood_ebb-forward/y-BE-'+str(BE)[:-2]
    df = h5py.File(output_dir+'/diagnostic_'+filenames[i]+'.hdf5','r+')
    # df = h5py.File(output_dir+'/diagnostic_sediment_errors.hdf5','r+')

    for name,data in df.items():
        if name == 'time':
            pass
        else:
            print(name)

            print(data[0],'\n',len(data),data[-1])
            if name == 'RMSE_current':
                a = int(len(data)/6)
                b = np.linspace(492,491+a,a)
                plt.plot(b,data[::6],'.-',label='4op-'+str(i))
                plt.legend()

BE = 0.0
for i in [1,2]:
    output_dir =  '../../../outputs/3.environment/discrete/flood_ebb-forward/y-BE-'+str(BE)[:-2]
    df = h5py.File(output_dir+'/diagnostic_'+filenames[i]+'.hdf5','r+')
    # df = h5py.File(output_dir+'/diagnostic_sediment_errors.hdf5','r+')

    for name,data in df.items():
        if name == 'time':
            pass
        else:
            print(name)

            print(data[0],'\n',len(data),data[-1])
            if name == 'RMSE_current':
                a = int(len(data)/6)
                b = np.linspace(492,491+a,a)
                plt.plot(b,data[::6],'.-',label='for-'+str(i))
                plt.legend()
plt.savefig('fig.jpg',dpi=300)

     
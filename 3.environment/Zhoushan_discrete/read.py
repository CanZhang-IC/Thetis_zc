import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
import os


thedir=os.getcwd()  

filenames = ['controls','velocity_errors','velocity_errors2','each_farm_optimisation','eachturbine','farm_optimisation','functional','turbine','volume2d']

BE = 8.0

output_dir =  '../../../outputs/2.economy/discrete/flood_ebb/cable-BE-'+str(BE)[:-2]
df = h5py.File(output_dir+'/diagnostic_'+filenames[-3]+'.hdf5','r+')
# df = h5py.File(output_dir+'/diagnostic_sediment_errors.hdf5','r+')

for name,data in df.items():
    if name == 'time':
        pass
    else:
        print(name)

        print(data[0],'\n',len(data),list(data[-1]))
        # if name == 'RMSE_current':
        a = int(len(data))
        b = np.linspace(492,491+a,a)
        plt.plot(b,data,'.-')
        plt.legend()

plt.savefig('fig.jpg',dpi=300)

     
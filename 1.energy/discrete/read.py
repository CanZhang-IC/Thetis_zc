import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
import os


thedir=os.getcwd()  

filenames = ['controls','derivatives','each_farm_optimisation','eachturbine','farm_optimisation','functional','turbine','volume2d']

output_dir =  '../../../outputs/1.energy/discrete/all-90'
df = h5py.File(output_dir+'/diagnostic_'+filenames[-3]+'.hdf5','r+')
# df = h5py.File(output_dir+'/diagnostic_sediment_errors.hdf5','r+')

for name,data in df.items():
    if name == 'time':
        pass
    else:
        print(name)
        # for i,ii in enumerate(data):
        #     print(i, ii)
        print(data[0],'\n',len(data),data[-1])
    a = len(data)
    b = np.linspace(1,a,a)
    plt.plot(b,data,'.-',)
plt.savefig('fig.jpg',dpi=300)

     
import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
import os


thedir=os.getcwd()  

filenames = ['controls','derivatives','each_farm_optimisation','eachturbine','farm_optimisation','functional','turbine','volume2d']

BE,BE_sediment = 5.0,25.0
output_dir =  '../../../outputs/3.environment/zhoushan-continuous-op/behind_notfrom0/'+str(BE)[:-2]+'_'+str(BE_sediment)[:-2]
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

     
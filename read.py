import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
import os


thedir=os.getcwd()  

filenames = ['controls','derivatives','each_farm_optimisation','eachturbine','farm_optimisation','functional','turbine','volume2d']

# df = h5py.File('../../outputs/ideal_cable/only_cable/diagnostic_'+filenames[-3]+'.hdf5','r+')
df = h5py.File('./diagnostic_'+filenames[-3]+'.hdf5','r+')
for name,data in df.items():
    for i,ii in enumerate(data):
        print(i, ii)
    # print('1',data[1])
    print(len(data)-1,data[-1])
    incremet = (float(data[-1])-float(data[0]))/float(data[0])*100
    print('The increment is {0:0.2f}%.'.format(incremet))
    a = len(data)
    b = np.linspace(1,a,a)
    plt.plot(b,data,'.-',)
plt.savefig('fig.jpg',dpi=300)

     
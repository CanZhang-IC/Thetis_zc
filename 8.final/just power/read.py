import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
import os


thedir=os.getcwd()  

filenames = ['controls','derivatives','each_farm_optimisation','eachturbine','farm_optimisation','functional','turbine','volume2d']

output_dir = '../../../outputs/8.final/just_power2'
df = h5py.File(output_dir+'/diagnostic_'+filenames[-3]+'.hdf5','r+')

for name,data in df.items():
    if name == 'time':
        pass
    else:
        print(name)
        print(list(data[0]))
        print(len(data),list(data[-1]))
        # for i,ii in enumerate(data):
        #     # if i ==0 :
        #     print(i, list(ii))
    # print('1',data[1])
#     print(len(data)-1,data[-1])
#     incremet = (float(data[-1])-float(data[0]))/float(data[0])*100
#     print('The increment is {0:0.2f}%.'.format(incremet))
        a = len(data)
        b = np.linspace(1,a,a)
        plt.plot(b,data,'.-',)
        # plt.ylim(0,2500)
    plt.savefig('fig.jpg',dpi=300)

     
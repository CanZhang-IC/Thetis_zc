import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../..')
import os
from prepare_cable.Hybrid_Code import CableCostGA

thedir=os.getcwd()  

filenames = ['controls','derivatives','each_farm_optimisation','eachturbine','farm_optimisation','functional','turbine','volume2d']

BE = 1.0

start_from_initial = True
output_dir = '../../../outputs/2.economy/discrete/intermediate-2-l_y/yaw-cable-BE-'+str(BE)[:-2]


# for iname in range(5,46,10):
df = h5py.File(output_dir+'/diagnostic_'+filenames[-4]+'.hdf5','r+')

for name,data in df.items():
    if name == 'average_power':
    #     pass
    # else:
        print(name)
        print(data[0])
        print(len(data),list(data[-1]))

        a = len(data)
        b = np.linspace(1,a,a)
        plt.plot(b,data,'.-',)
        # plt.ylim(0,2500)
    plt.savefig('fig.jpg',dpi=300)

df = h5py.File(output_dir+'/diagnostic_'+filenames[0]+'.hdf5','r+')

for name,data in df.items():
    all_controls = list(data[-1])
    turbine_locations = all_controls[:24]#[x for coord in turbine_location for x in coord]
    landpointlocation = [444000,3323000]
    CC = CableCostGA(turbine_locations, substation_location=landpointlocation)
    print ('Cable length:',CC.compute_cable_cost())



import datetime
from tkinter.font import families
from pylab import *
import matplotlib
import h5py
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from math import sqrt,atan
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score

zw = matplotlib.font_manager.FontProperties(fname='/media/can/can_disk/simsun.ttc',size=15)
# yw = matplotlib.font_manager.FontProperties(fname='/media/can/can_disk/times.ttf')
plt.rc('font',family='Times New Roman')


stationname='huluele'
print(stationname)

data_time=[]

start = datetime.datetime(2013,8,15,0,0,0)#start time
stop = datetime.datetime(2013,8,25,0,0,1)#end time
delta = datetime.timedelta(hours=1)
date = mpl.dates.drange(start,stop,delta)
data_time.append(date)

start = datetime.datetime(2013,8,15,0,0,0)#start time
stop = datetime.datetime(2013,8,25,0,0,1)#end time
delta = datetime.timedelta(hours=0.5)
date = mpl.dates.drange(start,stop,delta)
data_time.append(date)

start = datetime.datetime(2013,8,15,0,0,0)#start time
stop = datetime.datetime(2013,8,25,0,0,0)#end time
delta = datetime.timedelta(minutes=5)
date = mpl.dates.drange(start,stop,delta)
data_time.append(date)

fig, axs = plt.subplots(1, 1, figsize=(16,9))

### Read measured data ###
file_name = stationname + '.csv'
elevation = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(2),skiprows=1) 
# data formatting:
both = []
# remove empty values
for i in range(len(elevation)):
    if elevation[i] != '-':
        both.append([elevation[i]])
measured_elevation=[]
for i in range(len(both)):
    measured_elevation.append(float(both[i][0]))
axs.plot( data_time[0], measured_elevation[int(3*24):-int(2*24)], 'ko', label='实测数据')

thetisfilenames=[
        'continuous-4cores','continuous-4cores-tideforcing','discrete-4cores','discrete-4cores-tideforcing'
      ]
legend_names = ['$B_I$','$B_O$','$C_I$','$C_O$']
for ii,thetisfilename in enumerate(thetisfilenames):    
    ### Read Thetis data ###        
    det_file = '../../outputs/0.validation/'+thetisfilename+"/diagnostic_detectors.hdf5"
    df = h5py.File(det_file, 'r+')
    thetis_elevation=[]
    spin = int((3*24)*60*60/60/5)
    eight=int(8*60/5)
    for name, data in df.items():
        #print(name)
        if name == stationname:
            thetis_elevation.append(data[spin-eight:-eight,0])
    
    print(len(thetis_elevation[0]))
    # if thetisfilename == 'continuous-4cores-tideforcing':
    #     t_e_end = int(3*24*60/5)+2880
    # else:
    #     t_e_end = -int(2*24*60/5)
    t_e_end = -int(2*24*60/5)  
    axs.plot( data_time[2], thetis_elevation[0][int(3*24*60/5):t_e_end],  label=legend_names[ii])
    axs.set_xlabel('日期', fontsize=25,fontproperties=zw)
    axs.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    axs.set_ylabel('潮位 $[m]$', fontsize=25,fontproperties=zw)
    axs.legend(prop=zw)
    thetis_rmse_elevation1=mean_squared_error(np.array(thetis_elevation[0][int(3*24*60/5):t_e_end:12]),np.array(measured_elevation[int(3*24):-int(2*24)-1]), squared=False)
    thetis_r2_elevation1=r2_score(np.array(thetis_elevation[0][int(3*24*60/5):t_e_end:12]),np.array(measured_elevation[int(3*24):-int(2*24)-1]))
    print(thetisfilename)
    print(thetis_rmse_elevation1,thetis_r2_elevation1)

plt.xticks(data_time[0][::48],size=20)
plt.yticks(size=20)
# show()
plt.savefig('葫芦岛潮位BC.png', dpi =300)

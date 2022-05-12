import datetime
from pylab import *
import matplotlib
import h5py
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from math import sqrt,atan
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score


stationname='taohuaele'
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
axs.plot( data_time[0], measured_elevation[int(3*24):-int(2*24)], 'ko', label='Measurement')

# #read OTF data
# file_name = 'otf_water_level2.csv'
# otf_water_level = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(1),skiprows=0)  
# otf_x = []
# for i in range(len(otf_water_level)):
#     if otf_water_level[i] == None:
#         pass
#     else:
#         otf_x.append(float(otf_water_level[i])) 
# axs.plot(data_time[1], otf_x[int(24*60/30):], '-' ,  color='r',label='OTF_NF', linewidth=2)
# thetis_rmse_elevation2=mean_squared_error(np.array(otf_x[int(24*60/30)::2]),np.array(measured_elevation[int(3*24):-int(2*24)]), squared=False)
# thetis_r2_elevation2=r2_score(np.array(otf_x[int(24*60/30)::2]),np.array(measured_elevation[int(3*24):-int(2*24)]))
# print('OTF_near-field')
# print(thetis_rmse_elevation2,thetis_r2_elevation2)

# ### Far-field B5 Read OTF data ###
# file_name = 'otf_water_level.csv'
# otf_water_level = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(0),skiprows=1)  
# otf_x = []
# for i in range(len(otf_water_level)):
#     otf_x.append(float(otf_water_level[i])) 
# axs.plot(data_time[1], otf_x[80:-16], '-' , color='orange',label='OTF_FF', linewidth=2)
# thetis_rmse_elevation2=mean_squared_error(np.array(otf_x[80:-16:2]),np.array(measured_elevation[int(3*24):-int(2*24)]), squared=False)
# thetis_r2_elevation2=r2_score(np.array(otf_x[80:-16:2]),np.array(measured_elevation[int(3*24):-int(2*24)]))
# print('OTF_far-field')
# print(thetis_rmse_elevation2,thetis_r2_elevation2)

### Far-field A4 Read OTF data ###
file_name = 'otf_water_level.csv'
otf_water_level = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(1),skiprows=1)  
otf_x = []
for i in range(len(otf_water_level)):
    otf_x.append(float(otf_water_level[i])) 
axs.plot(data_time[1], otf_x[80:-16], '-' , color='orange',label='OTF_FF', linewidth=2)
thetis_rmse_elevation2=mean_squared_error(np.array(otf_x[80:-16:2]),np.array(measured_elevation[int(3*24):-int(2*24)]), squared=False)
thetis_r2_elevation2=r2_score(np.array(otf_x[80:-16:2]),np.array(measured_elevation[int(3*24):-int(2*24)]))
print('OTF_far-field')
print(thetis_rmse_elevation2,thetis_r2_elevation2)



thetisfilenames=['onemin-8cores-huluthreepart']

for thetisfilename in thetisfilenames:    
    ### Read Thetis data ###        
    det_file = "../outputs/"+thetisfilename+"/diagnostic_detectors.hdf5"
    df = h5py.File(det_file, 'r+')
    thetis_elevation=[]
    spin = int((3*24)*60*60/60/5)
    eight=int(8*60/5)
    for name, data in df.items():
        #print(name)
        if name == stationname:
            thetis_elevation.append(data[spin-eight:-eight,0])

    axs.plot( data_time[2], thetis_elevation[0][int(3*24*60/5):-int(2*24*60/5)],'b-',  label='Thetis-5min')
    axs.set_xlabel('$Time$ $[date]$', fontsize=25)
    axs.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    axs.set_xticks(data_time[0][::48])
    axs.set_ylabel('$\eta$ $[m]$', fontsize=25)
    axs.xaxis.set_tick_params(labelsize=20)
    axs.yaxis.set_tick_params(labelsize=20)
    axs.legend(loc='best', fontsize=15)
thetis_rmse_elevation1=mean_squared_error(np.array(thetis_elevation[0][int(3*24*60/5):-int(2*24*60/5):12]),np.array(measured_elevation[int(3*24):-int(2*24)-1]), squared=False)
thetis_r2_elevation1=r2_score(np.array(thetis_elevation[0][int(3*24*60/5):-int(2*24*60/5):12]),np.array(measured_elevation[int(3*24):-int(2*24)-1]))
print('Thetis_5min')
print(thetis_rmse_elevation1,thetis_r2_elevation1)


thetisfilenames=['vsotf']

for thetisfilename in thetisfilenames:    
    ### Read Thetis data ###        
    det_file = "../outputs/"+thetisfilename+"/diagnostic_detectors.hdf5"
    df = h5py.File(det_file, 'r+')
    thetis_elevation=[]
    for name, data in df.items():
        #print(name)
        if name == stationname:
            thetis_elevation.append(data[:,0])

    axs.plot( data_time[1], thetis_elevation[0][79:-16],'g-',  label='Thetis-30min')
    axs.set_xlabel('$Time$ $[date]$', fontsize=25)
    axs.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    axs.set_xticks(data_time[0][::48])
    axs.set_ylabel('$\eta$ $[m]$', fontsize=25)
    axs.xaxis.set_tick_params(labelsize=20)
    axs.yaxis.set_tick_params(labelsize=20)
    axs.legend(loc='best', fontsize=15)
thetis_rmse_elevation1=mean_squared_error(np.array(thetis_elevation[0][79:-16:2]),np.array(measured_elevation[int(3*24):-int(2*24)]), squared=False)
thetis_r2_elevation1=r2_score(np.array(thetis_elevation[0][79:-16:2]),np.array(measured_elevation[int(3*24):-int(2*24)]))
print('Thetis_30min')
print(thetis_rmse_elevation1,thetis_r2_elevation1)
show()
#plt.savefig('taohua-elevation.png', dpi =300)

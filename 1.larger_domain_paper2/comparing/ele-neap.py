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

data_time=[]

start = datetime.datetime(2013,8,16,10,0,0)
stop = datetime.datetime(2013,8,17,11,0,1)
delta = datetime.timedelta(hours=1)
date = mpl.dates.drange(start,stop,delta)
data_time.append(date)

###2013/8/15  me:3*24 oe:24*60/30 te:3*24*60/5
###2013/8/25  me:-2*24 oe: te:-2*24*60/5
time_begin = 1*24+10
time_end   = 7*24+13

fig, axs = plt.subplots(1, 1, figsize=(16,9))

### Read measured data ###
file_name = 'taohuaele.csv'
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

#read OTF data
file_name = 'otf_water_level2.csv'
otf_water_level = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(1),skiprows=0) 
     
otf_x = []
for i in range(len(otf_water_level)):
    if otf_water_level[i] == None:
        pass
    else:
        otf_x.append(float(otf_water_level[i])) 





### Read Thetis data ###  
thetisfilenames=['onemin-8cores-huluthreepart']
for thetisfilename in thetisfilenames:    
          
    det_file = "../outputs/"+thetisfilename+"/diagnostic_detectors.hdf5"
    df = h5py.File(det_file, 'r+')
    thetis_elevation=[]
    spin = int((3*24)*60*60/60/5)
    eight=int(8*60/5)
    for name, data in df.items():
        #print(name)
        if name == 'taohuaele':
            thetis_elevation.append(data[spin-eight:-eight,0])
    #axs.plot( data_time[2], thetis_elevation[0][int(3*24*60/5):-int(2*24*60/5)],'.-',  label='Thetis-thetis_elevation')
   

me=measured_elevation[int(3*24+time_begin):-int(2*24+time_end)]
oe=otf_x[int((24+time_begin)*60/30):-int(time_end*60/30):2]
otf_rmse_elevation=mean_squared_error(np.array(oe),np.array(me), squared=False)
otf_r2_elevation=r2_score(np.array(oe),np.array(me))
thetis_min_rmse=10
thetis_min_r2=10
a_min=0
j_min=0
for a in [12,6,4,3,2]:
    for j in range(a):
        te1=[]
        for i in range(int((3*24+time_begin)*60/5),len(thetis_elevation[0])-int((2*24+time_end)*60/5)+2,a):
            te1.append(np.mean(thetis_elevation[0][i-j:i+a-j]))
        rmse=mean_squared_error(np.array(te1[::int(12/a)]),np.array(me), squared=False)
        r2  =r2_score(np.array(te1[::int(12/a)]),np.array(me))
        if rmse < thetis_min_rmse:
            thetis_min_rmse = rmse 
            thetis_min_r2   = r2
            a_min=a
            j_min=j

print('otf rmse:',otf_rmse_elevation,'thetis rmse:',thetis_min_rmse)
print('otf r2:',otf_r2_elevation,'thetis r2:',thetis_min_r2)
print(a_min,j_min)

axs.plot( data_time[0], me, 'ko', label='measurement-A2')
axs.plot(data_time[0], oe, 'r.-', label='OpenTidalFarm-A2', linewidth=2)

te=[]
for i in range(int((3*24+time_begin)*60/5),len(thetis_elevation[0])-int((2*24+time_end)*60/5)+2,a_min):
    te.append(np.mean(thetis_elevation[0][i-j_min:i+a_min-j_min]))
axs.plot(data_time[0], te[::int(12/a_min)], 'b.-', label='Thetis-A2', linewidth=2)

axs.set_xlabel('$time$', fontsize=14)
axs.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
axs.set_xticks(data_time[0][::12])
axs.set_ylabel('$elevation(m)$', fontsize=14)
axs.legend(loc='best', fontsize=9)
#show()
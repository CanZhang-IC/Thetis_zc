import datetime
from pylab import *
import matplotlib
import h5py
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from math import sqrt,atan

#read measured and otf data
file_name = 'A5.csv'
velocity = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(1),skiprows=1) # velocity in 2nd column
otf_velocity = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(3),skiprows=1)
# data formatting:
both = []
# remove empty values
for i in range(len(time)):
    if velocity[i] != '-':
        both.append([velocity[i],otf_velocity[i]])
x_neap = []
otf_x_neap = []
for i in range(51):
    x_neap.append(float(both[i][0])) # converting water level from cm to m
    otf_x_neap.append(float(both[i][1]))
x_meso = []
otf_x_meso = []
for i in range(51,102):
    x_meso.append(float(both[i][0])) # converting water level from cm to m
    otf_x_meso.append(float(both[i][1]))
x_spring = []
otf_x_spring = []
for i in range(102,153):
    x_spring.append(float(both[i][0])) # converting water level from cm to m
    otf_x_spring.append(float(both[i][1]))


#read Thetis data
det_file = "../outputs/v1-m0.02-10/diagnostic_detectors.hdf5"
df = h5py.File(det_file, 'r+')

xvelocity=[]
yvelocity=[]
for name, data in df.items():
    if name == 'A5':       
        xvelocity.append(data[:,1])
        yvelocity.append(data[:,2])

thetis_velocity=[]
thetis_veldirection=[]
for i in range(len(xvelocity[0])):
    thetis_velocity.append(sqrt(xvelocity[0][i]**2+yvelocity[0][i]**2))
    if yvelocity[0][i] < 0 :
        thetis_veldirection.append(atan(xvelocity[0][i]/yvelocity[0][i])/np.pi*180+180)
    else:
        if xvelocity[0][i] < 0:
            thetis_veldirection.append(atan(xvelocity[0][i]/yvelocity[0][i])/np.pi*180+360)
        else:
            thetis_veldirection.append(atan(xvelocity[0][i]/yvelocity[0][i])/np.pi*180)

thetis_x_all = []
for i in range(len(thetis_velocity)):
    '''
    if thetis_veldirection[i]> 200:
        thetis_x_all.append(-float(thetis_velocity[i]))
    else:
        thetis_x_all.append(float(thetis_velocity[i]))
    '''
    thetis_x_all.append(float(thetis_velocity[i]))



startneap = datetime.datetime(2013,8,9,8,0,0)#起始时间
stopneap = datetime.datetime(2013,8,27,8,0,0)#停止时间
delta = datetime.timedelta(seconds=60)
dateall = mpl.dates.drange(startneap,stopneap,delta)

startneap = datetime.datetime(2013,8,16,10,0,0)#起始时间
stopneap = datetime.datetime(2013,8,17,11,0,0)#停止时间
delta = datetime.timedelta(minutes=10)
dateneap = mpl.dates.drange(startneap,stopneap,delta)

startmeso = datetime.datetime(2013,8,19,14,0,0)#起始时间
stopmeso = datetime.datetime(2013,8,20,15,0,1)#停止时间
delta = datetime.timedelta(minutes=10)
datemeso = mpl.dates.drange(startmeso,stopmeso,delta)

startspring = datetime.datetime(2013,8,23,10,0,0)#起始时间
stopspring = datetime.datetime(2013,8,24,11,0,1)#停止时间
delta = datetime.timedelta(minutes=10)
dataspring = mpl.dates.drange(startspring,stopspring,delta)



'''
fig, axs = plt.subplots(1, 1, figsize=(10, 10))
axs.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
axs.set_xticks(dateall[0:-1:5*24*60])
axs.plot(dateall, thetis_x_all,label='finemesh')
axs.legend(loc='best', fontsize=14)
axs.plot( dateneap,x_neap, 'ko')
axs.plot( dateneap,otf_x_neap, 'g')
axs.plot(datemeso, x_meso, 'ko')
axs.plot( dataspring, x_spring,'ko')
axs.plot(datemeso, otf_x_meso, 'g')
axs.plot( dataspring, otf_x_spring,'g')
'''



# set up figure with 4 subplots to plot 4 resolutions
fig, axs = plt.subplots(3, 1, figsize=(10, 10))
# reshape so that we can iterate below over axs[i] instead of ax[i,j]
axs = axs.reshape(-1)
fig.tight_layout(h_pad=4)
zz=zeros(len(dateneap))


axs[0].plot( dateneap,x_neap, 'ko', label='measurement')
axs[0].plot( dateall[10200:11701],thetis_x_all[10200:11701], 'b', label='Thetis')
axs[0].plot( dateneap,otf_x_neap, 'g', label='OpenTidalFarm')
axs[0].plot( dateneap,zz, 'r--')
axs[0].set_xlabel('$time(s)$', fontsize=14)
axs[0].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
axs[0].set_xticks(dateneap[0:-1:12])
axs[0].set_ylabel('$velocity(m/s)$', fontsize=14)
axs[0].set_title('Comparisons between simulated and measured velocity', fontsize=14)
axs[0].legend(loc='best', fontsize=14)


axs[1].plot(datemeso, x_meso, 'ko', label='measurement')
axs[1].plot(dateall[14760:16261],thetis_x_all[14760:16261], 'b', label='Thetis')
axs[1].plot(datemeso, zz, 'r--')
axs[1].plot(datemeso, otf_x_meso, 'g', label='OpenTidalFarm')
axs[1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
axs[1].set_xticks(datemeso[0:-1:12])
axs[1].set_xlabel('$time(s)$', fontsize=14)
axs[1].set_ylabel('$velocity(m/s)$', fontsize=14)
axs[1].set_title('Comparisons between simulated and measured velocity', fontsize=14)
axs[1].legend(loc='best', fontsize=14)


axs[2].plot( dataspring, x_spring,'ko', label='measurement')
axs[2].plot(dateall[20280:21781],thetis_x_all[20280:21781], 'b', label='Thetis')
axs[2].plot(dataspring, otf_x_spring, 'g', label='OpenTidalFarm')
axs[2].plot(dataspring, zz, 'r--')
axs[2].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
axs[2].set_xticks(dataspring[0:-1:12])
axs[2].set_xlabel('$time(s)$', fontsize=14)
axs[2].set_ylabel('$velocity(m/s)$', fontsize=14)
axs[2].set_title('Comparisons between simulated and measured velocity', fontsize=14)
axs[2].legend(loc='best', fontsize=14)




show()

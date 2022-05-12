from uptide import *
import datetime
from pandas.plotting import register_matplotlib_converters
import math
from pylab import *
import pandas as pd

import matplotlib
import h5py
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from math import sqrt,atan
import shutil
#read measured data
file_name = '5_velocity_measured_otf.csv'
direction = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(3),skiprows=1) # direction in 2nd column
otf_direction = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(5),skiprows=1)
# data formatting:
both = []
# remove empty values
for i in range(len(direction)):
    if direction[i] != '-':
        both.append([direction[i],otf_direction[i]])
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

fig, axs = plt.subplots(3, 1, figsize=(10, 10))
# reshape so that we can iterate below over axs[i] instead of ax[i,j]
axs = axs.reshape(-1)
fig.tight_layout(h_pad=4)
startneap = datetime.datetime(2013,8,16,10,0,0)#起始时间
stopneap = datetime.datetime(2013,8,17,11,0,0)#停止时间
delta = datetime.timedelta(hours=0.5)
t_neap = mpl.dates.drange(startneap,stopneap,delta)

startmeso = datetime.datetime(2013,8,19,14,0,0)#起始时间
stopmeso = datetime.datetime(2013,8,20,15,0,1)#停止时间
delta = datetime.timedelta(hours=0.5)
t_meso = mpl.dates.drange(startmeso,stopmeso,delta)

startspring = datetime.datetime(2013,8,23,10,0,0)#起始时间
stopspring = datetime.datetime(2013,8,24,11,0,1)#停止时间
delta = datetime.timedelta(hours=0.5)
t_spring = mpl.dates.drange(startspring,stopspring,delta)

axs[0].plot( t_neap,x_neap, 'ko', label='measurement')
axs[0].plot( t_neap,otf_x_neap, 'r', label='OpenTidalFarm',linewidth=3)

axs[1].plot(t_meso, x_meso, 'ko', label='measurement')
axs[1].plot(t_meso, otf_x_meso, 'r', label='OpenTidalFarm',linewidth=3)

axs[2].plot( t_spring, x_spring,'ko', label='measurement')
axs[2].plot(t_spring, otf_x_spring, 'r', label='OpenTidalFarm',linewidth=3)

startneap = datetime.datetime(2013,8,16,10,0,0)#起始时间
stopneap = datetime.datetime(2013,8,17,11,0,0)#停止时间
delta = datetime.timedelta(seconds=60)
dateneap = mpl.dates.drange(startneap,stopneap,delta)

startmeso = datetime.datetime(2013,8,19,14,0,0)#起始时间
stopmeso = datetime.datetime(2013,8,20,15,0,1)#停止时间
delta = datetime.timedelta(seconds=60)
datemeso = mpl.dates.drange(startmeso,stopmeso,delta)

startspring = datetime.datetime(2013,8,23,10,0,0)#起始时间
stopspring = datetime.datetime(2013,8,24,11,0,1)#停止时间
delta = datetime.timedelta(seconds=60)
dataspring = mpl.dates.drange(startspring,stopspring,delta)

names=['onemin-8cores-huluthreepart']

for filename in names:
    det_file = "../outputs/"+thetisfilename+"/diagnostic_detectors.hdf5"
    df = h5py.File(det_file, 'r+')
    xvelocity=[]
    yvelocity=[]
    for name, data in df.items():
        if name == 'B5':
            print(name)
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
    thetis_t_all = np.arange(0,1555200,60)
    for i in range(len(thetis_velocity)):
            thetis_x_all.append(float(thetis_veldirection[i]))

    thetis_x_neap = []
    #thetis_t_neap = np.arange(381600,471600,1800)
    for i in range(1501):
        thetis_x_neap.append(thetis_x_all[10200+i]) # converting water level from cm to m

    thetis_x_meso = []
    #thetis_t_meso = np.arange(655200,745200,1800)
    for i in range(1500,3001):
        thetis_x_meso.append(thetis_x_all[14760+i]) # converting water level from cm to m

    thetis_x_spring = []
    #thetis_t_spring = np.arange(986400,1076400,1800)
    for i in range(3000,4501):
        thetis_x_spring.append(thetis_x_all[20280+i]) # converting water level from cm to m



    axs[0].plot( dateneap,thetis_x_neap, label=filename)
    axs[0].set_xlabel('$time(s)$', fontsize=14)
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
    axs[0].set_xticks(dateneap[0:-1:12*30])
    axs[0].set_ylabel('$direction(°)$', fontsize=14)
    axs[0].set_title('Comparisons between simulated and measured direction', fontsize=14)
    axs[0].legend(loc='best', fontsize=14)



    axs[1].plot(datemeso, thetis_x_meso, label=filename)
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
    axs[1].set_xticks(datemeso[0:-1:12*30])
    axs[1].set_xlabel('$time(s)$', fontsize=14)
    axs[1].set_ylabel('$direction(°)$', fontsize=14)
    axs[1].set_title('Comparisons between simulated and measured direction', fontsize=14)
    axs[1].legend(loc='best', fontsize=6)



    axs[2].plot(dataspring, thetis_x_spring, label=filename)
    axs[2].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
    axs[2].set_xticks(dataspring[0:-1:12*30])
    axs[2].set_xlabel('$time(s)$', fontsize=14)
    axs[2].set_ylabel('$direction(°)$', fontsize=14)
    axs[2].set_title('Comparisons between simulated and measured direction', fontsize=14)
    axs[2].legend(loc='best', fontsize=6)




show()


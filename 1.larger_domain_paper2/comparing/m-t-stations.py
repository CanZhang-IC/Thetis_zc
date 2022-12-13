import datetime
from pylab import *
import matplotlib
import h5py
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from math import sqrt,atan

stationnames=['A4','B5']#'A1','A2','A3','A4','A5','B1','B2','B3','B4','B5']
for stationname in stationnames:
    print(stationname)
    ### Read measured data ###
    file_name = stationname+'.csv'
    velocity = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(1),skiprows=1) 
    # data formatting:
    measured_velocity = []
    # remove empty values
    for i in range(len(velocity)):
        if velocity[i] != '-':
            measured_velocity.append(float(velocity[i]))

    ### Read Thetis data ###
    name='minimumdepth(2)-manning(0.02-0.3)-viscosity(5-1000)-simulation-test'
    det_file = "../outputs/"+name+"/diagnostic_detectors.hdf5"
    df = h5py.File(det_file, 'r+')
    xvelocity=[]
    yvelocity=[]
    for name, data in df.items():
        if name == stationname:
        	xvelocity.append(data[:,1])
            yvelocity.append(data[:,2])
    thetis_velocity=[]
    for i in range(len(xvelocity[0])):
    	thetis_velocity.append(sqrt(xvelocity[0][i]**2+yvelocity[0][i]**2))
    
    thetis_velocity3=[]
    for i in range(2040,2341,2):
        thetis_velocity3.append(np.mean(thetis_velocity[i:i+2]))
    for i in range(2952,3253,2):
        thetis_velocity3.append(np.mean(thetis_velocity[i:i+2]))
    for i in range(4056,4357,2):
        thetis_velocity3.append(np.mean(thetis_velocity[i:i+2]))
    

    #x-time('%Y-%m-%d %H:%M:%S')
    data_time=[]

    startneap = datetime.datetime(2013,8,16,10,0,0)#start time
    stopneap = datetime.datetime(2013,8,17,11,0,0)#end time
    delta = datetime.timedelta(minutes=10)
    dateneap = mpl.dates.drange(startneap,stopneap,delta)
    data_time.append(dateneap)

    startmeso = datetime.datetime(2013,8,19,14,0,0)
    stopmeso = datetime.datetime(2013,8,20,15,0,1)
    delta = datetime.timedelta(minutes=10)
    datemeso = mpl.dates.drange(startmeso,stopmeso,delta)
    data_time.append(datemeso)

    startspring = datetime.datetime(2013,8,23,10,0,0)
    stopspring = datetime.datetime(2013,8,24,11,0,1)
    delta = datetime.timedelta(minutes=10)
    dataspring = mpl.dates.drange(startspring,stopspring,delta)
    data_time.append(dataspring)

    fig, axs = plt.subplots(3, 1, figsize=(10, 10))
    axs = axs.reshape(-1) # reshape so that we can iterate below over axs[i] instead of ax[i,j]
    fig.tight_layout(pad=3, h_pad=1.5)
    for i in range(3):
        axs[i].plot( data_time[i], measured_velocity[151*i:151*(i+1)], 'ko', label='measurement-'+stationname)
        axs[i].plot( data_time[i], thetis_velocity3[151*i:151*(i+1)], 'b', label='Thetis-'+stationname)
        axs[i].set_xlabel('$time$', fontsize=14)
        axs[i].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
        axs[i].set_xticks(data_time[i][0:-1:24])
        axs[i].set_ylabel('$velocity(m/s)$', fontsize=14)
        axs[i].legend(loc='lower right', fontsize=9)

show()

import datetime
from pylab import *
import matplotlib
import h5py
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from math import sqrt,atan



stationnames=['B1','B2','B3','B4','B5']
stationnames=['B5']
for stationname in stationnames:
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
    fig, axs = plt.subplots(3, 1, figsize=(16,9))
    axs = axs.reshape(-1) # reshape so that we can iterate below over axs[i] instead of ax[i,j]
    fig.tight_layout(pad=4,h_pad=3.5)

    otf_time = []

    startneap = datetime.datetime(2013,8,16,10,0,0)#start time
    stopneap = datetime.datetime(2013,8,17,11,0,0)#end time
    delta = datetime.timedelta(minutes=30)
    dateneap = mpl.dates.drange(startneap,stopneap,delta)
    otf_time.append(dateneap)

    startmeso = datetime.datetime(2013,8,19,14,0,0)
    stopmeso = datetime.datetime(2013,8,20,15,0,1)
    delta = datetime.timedelta(minutes=30)
    datemeso = mpl.dates.drange(startmeso,stopmeso,delta)
    otf_time.append(datemeso)

    startspring = datetime.datetime(2013,8,23,10,0,0)
    stopspring = datetime.datetime(2013,8,24,11,0,1)
    delta = datetime.timedelta(minutes=30)
    dataspring = mpl.dates.drange(startspring,stopspring,delta)
    otf_time.append(dataspring)

    ### Read measured data ###
    file_name =stationname+'.csv'
    velocity = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(1),skiprows=1) # velocity in 2nd column
    direction = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(2),skiprows=1)

    # data formatting:
    both = []
    # remove empty values
    for i in range(len(velocity)):
        if velocity[i] != '-':
            both.append([velocity[i], direction[i]])

    velocity_magnitude = []
    velocity_direction= []

    for i in range(len(both)):
        velocity_magnitude.append(float(both[i][0]))
        velocity_direction.append(float(both[i][1]))
    
    for i in range(3):
        axs[i].plot( data_time[i], velocity_direction[151*i:151*(i+1)], 'ko', label='Measurement')


        
    # ### Read OTF data ###
    # file_name = 'v_otf.csv'
    # otf_v = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(3),skiprows=0) 
    # # data formatting:
    # all_otf_v = []
    # # remove empty values
    # for i in range(len(otf_v)):
    #     if otf_v[i] != '-':
    #         all_otf_v.append(float(otf_v[i]))

    # axs[0].plot(otf_time[0], all_otf_v[149:200], '.-' , color='r',label='OTF_NF')
    # axs[1].plot(otf_time[1], all_otf_v[300:351], '.-' , color='r', label='OTF_NF')
    # axs[2].plot(otf_time[2], all_otf_v[484:535], '.-' , color='r', label='OTF_NF')

    # ### Far-field B5 Read OTF data ###
    # file_name = 'v_otf.csv'
    # otf_xv = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(5),skiprows=0) 
    # otf_yv = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(6),skiprows=0)
    # # data formatting:
    # all_otf_v = []
    # # remove empty values
    # for i in range(len(otf_xv)):
    #     one_direction = np.arctan2(float(otf_xv[i]),float(otf_yv[i]))* 180 / np.pi
    #     if one_direction < 0:
    #         all_otf_v.append(one_direction+360)
    #     else:
    #         all_otf_v.append(one_direction)

    # axs[0].plot(otf_time[0], all_otf_v[147:198], '.-' , color='orange',label='OTF_FF')
    # axs[1].plot(otf_time[1], all_otf_v[298:349], '.-' , color='orange', label='OTF_FF')
    # axs[2].plot(otf_time[2], all_otf_v[482:533], '.-' , color='orange', label='OTF_FF')

    # ### Far-field A4 Read OTF data ###
    # file_name = 'v_otf.csv'
    # otf_xv = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(8),skiprows=0) 
    # otf_yv = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(9),skiprows=0)
    # # data formatting:
    # all_otf_v = []
    # # remove empty values
    # for i in range(len(otf_xv)):
    #     one_direction = np.arctan2(float(otf_xv[i]),float(otf_yv[i]))* 180 / np.pi
    #     if one_direction < 0:
    #         all_otf_v.append(one_direction+360)
    #     elif one_direction < 140:
    #         all_otf_v.append(360-one_direction)
    #     else:
    #         all_otf_v.append(one_direction)

    # axs[0].plot(otf_time[0], all_otf_v[147:198], '.-' , color='orange',label='OTF_FF')
    # axs[1].plot(otf_time[1], all_otf_v[298:349], '.-' , color='orange', label='OTF_FF')
    # axs[2].plot(otf_time[2], all_otf_v[482:533], '.-' , color='orange', label='OTF_FF')
        

    thetisfilenames=[
        'paper2validation','redata_5min_plus1',
        'redata_5min_normaldepth'
      ]
   

    for (ii,thetisfilename) in enumerate(thetisfilenames):    
        ### Read Thetis data ###        
        det_file = "../outputs/"+thetisfilename+"/diagnostic_detectors.hdf5"
        df = h5py.File(det_file, 'r+')
        xvelocity=[]
        yvelocity=[]
        for name, data in df.items():
            if name == stationname:       
                xvelocity.append(data[:,1])
                yvelocity.append(data[:,2])

        thetis_xvelocity3=[]
        thetis_yvelocity3=[]
        if thetisfilename == 'redata-dgdg-normaldepth' or thetisfilename == 'redata-dgdg-Latitude-30':
            for i in range(340,392):
                thetis_xvelocity3.append(np.mean(xvelocity[0][i]))
                thetis_yvelocity3.append(np.mean(yvelocity[0][i]))
            for i in range(492,544):
                thetis_xvelocity3.append(np.mean(xvelocity[0][i]))
                thetis_yvelocity3.append(np.mean(yvelocity[0][i]))
            for i in range(677,729):
                thetis_xvelocity3.append(np.mean(xvelocity[0][i]))
                thetis_yvelocity3.append(np.mean(yvelocity[0][i]))
        else:
            for i in range(2040,2341,2):
                thetis_xvelocity3.append(np.mean(xvelocity[0][i:i+2]))
                thetis_yvelocity3.append(np.mean(yvelocity[0][i:i+2]))
            for i in range(2952,3253,2):
                thetis_xvelocity3.append(np.mean(xvelocity[0][i:i+2]))
                thetis_yvelocity3.append(np.mean(yvelocity[0][i:i+2]))
            for i in range(4056,4357,2):
                thetis_xvelocity3.append(np.mean(xvelocity[0][i:i+2]))
                thetis_yvelocity3.append(np.mean(yvelocity[0][i:i+2]))
            
        thetis_velocity=[]
        thetis_direction=[]
        for i in range(len(thetis_xvelocity3)):
            thetis_velocity.append(sqrt(thetis_xvelocity3[i]**2+thetis_yvelocity3[i]**2))
            one_direction = np.arctan2(thetis_xvelocity3[i],thetis_yvelocity3[i])* 180 / np.pi
            if one_direction < 0:
                thetis_direction.append(one_direction + 360)
            else:
                thetis_direction.append(one_direction)


    
        
        # for i in range(len(thetis_direction)):
        #     if thetis_direction[i] < 140:
        #         thetis_direction[i] = 360 - thetis_direction[i]
        #     if velocity_direction[i] < 140:
        #         velocity_direction[i] = 360 - velocity_direction[i]
        
        tidenames = ['Neap','Intermediate','Spring']
        for i in range(3):
            if thetisfilename == 'redata-dgdg-normaldepth' or thetisfilename == 'redata-dgdg-Latitude-30':
                axs[i].plot( otf_time[i], thetis_direction[51*i:51*(i+1)],  label=thetisfilename)
            else:
                axs[i].plot(otf_time[i], thetis_direction[151*i:151*(i+1):3], '.-' ,label=thetisfilename)
            
            axs[2].set_xlabel('$Time$ $[date]$', fontsize=18)
            axs[i].set_ylabel('$\\theta$ $[°]$', fontsize=18)
            axs[i].xaxis.set_tick_params(labelsize=15)
            axs[i].yaxis.set_tick_params(labelsize=15)
            axs[i].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
            axs[i].set_xticks(data_time[i][0:-1:48])
            yy=np.linspace(0,360,7)
            axs[i].set_yticks(yy)
            #axs[i].set_ylim((0,360))
            axs[i].set_title('Velocity direction time series at B2 -- '+tidenames[i]+' tide',fontsize=18)
            axs[2].legend(loc='lower right', fontsize=14.5,ncol=4)
        #plt.savefig('taohua-direction.png', dpi =300)
show()


import datetime
from pylab import *
import matplotlib
import h5py
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from math import sqrt,atan



#stationnames=['A1','A2','A3','A4','A5','B1','B2','B3','B4','B5']
stationnames=['B5',]#'B5']#'A1','A2','A3','A4','A5','B1','B2','B3','B4','B5']
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

    startspring = datetime.datetime(2013,8,9,8,0,0)
    stopspring = datetime.datetime(2013,8,27,8,0,0)
    delta = datetime.timedelta(minutes=5)
    dataspring = mpl.dates.drange(startspring,stopspring,delta)
    data_time.append(dataspring)

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



    fig, axs = plt.subplots(3, 1, figsize=(16,9))
    axs = axs.reshape(-1) # reshape so that we can iterate below over axs[i] instead of ax[i,j]
    fig.tight_layout(pad=4,h_pad=3.5)

    ### Read measured data ###
    file_name = stationname+'.csv'
    velocity = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(1),skiprows=1) 
    # data formatting:
    measured_velocity = []
    # remove empty values
    for i in range(len(velocity)):
        if velocity[i] != '-':
            measured_velocity.append(float(velocity[i]))
    
    for i in range(3):
        axs[i].plot( data_time[i], measured_velocity[151*i:151*(i+1)], 'ko', label='Measurement')


    

    ### Read OTF data ###
    file_name = 'v_otf.csv'
    otf_v = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(2),skiprows=0) 
    # data formatting:
    all_otf_v = []
    # remove empty values
    for i in range(len(otf_v)):
        if otf_v[i] != '-':
            all_otf_v.append(float(otf_v[i]))

    # axs[0].plot(otf_time[0], all_otf_v[149:200], 'r', label='OTF_NF')
    # axs[1].plot(otf_time[1], all_otf_v[300:351], 'r',label='OTF_NF')
    # axs[2].plot(otf_time[2], all_otf_v[484:535], 'r',label='OTF_NF')

    # ### Far-field B5 Read OTF data ###
    # file_name = 'v_otf.csv'
    # otf_vx = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(5),skiprows=0) 
    # otf_vy = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(6),skiprows=0) 
    # # data formatting:
    # all_otf_v = []
    # # remove empty values
    # for i in range(len(otf_vx)):
    #     all_otf_v.append(sqrt(float(otf_vx[i])**2+float(otf_vy[i])**2))     
    # axs[0].plot(otf_time[0], all_otf_v[147:198], 'orange', label='OTF_FF')
    # axs[1].plot(otf_time[1], all_otf_v[298:349], 'orange',label='OTF_FF')
    # axs[2].plot(otf_time[2], all_otf_v[482:533], 'orange',label='OTF_FF')

    # ### Far-field A4 Read OTF data ###
    # file_name = 'v_otf.csv'
    # otf_vx = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(8),skiprows=0) 
    # otf_vy = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(9),skiprows=0) 
    # # data formatting:
    # all_otf_v = []
    # # remove empty values
    # for i in range(len(otf_vx)):
    #     all_otf_v.append(sqrt(float(otf_vx[i])**2+float(otf_vy[i])**2))     
    # axs[0].plot(otf_time[0], all_otf_v[147:198], 'orange', label='OTF_FF')
    # axs[1].plot(otf_time[1], all_otf_v[298:349], 'orange',label='OTF_FF')
    # axs[2].plot(otf_time[2], all_otf_v[482:533], 'orange',label='OTF_FF')
       
        

    thetisfilenames=[
        '5min-1cores-sediment'
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
        thetis_velocity=[]
        
        for i in range(len(xvelocity[0])):
            thetis_velocity.append(sqrt(xvelocity[0][i]**2+yvelocity[0][i]**2))
        

        
        thetis_velocity3=[]
        
        if thetisfilename=='onemin-8cores':
            for i in range(10200,11701,10):
                thetis_velocity3.append(np.mean(thetis_velocity[i:i+10]))
            for i in range(14760,16261,10):
                thetis_velocity3.append(np.mean(thetis_velocity[i:i+10]))
            for i in range(20280,21781,10):
                thetis_velocity3.append(np.mean(thetis_velocity[i:i+10]))
        elif thetisfilename == 'vsotf_otfmesh':
            for i in range(147,198):
                thetis_velocity3.append(thetis_velocity[i])
            for i in range(298,349):
                thetis_velocity3.append(thetis_velocity[i])
            for i in range(482,533):
                thetis_velocity3.append(thetis_velocity[i])


        else:
            for i in range(2040,2341,2):
                thetis_velocity3.append(np.mean(thetis_velocity[i:i+2]))
            for i in range(2952,3253,2):
                thetis_velocity3.append(np.mean(thetis_velocity[i:i+2]))
            for i in range(4056,4357,2):
                thetis_velocity3.append(np.mean(thetis_velocity[i:i+2]))     
        
        # for i in range(2040,2341,2):
            
        #     thetis_velocity3.append(thetis_velocity[i])
        # for i in range(2952,3253,2):
            
        #     thetis_velocity3.append(thetis_velocity[i])
        # for i in range(4056,4357,2):
            
        #     thetis_velocity3.append(thetis_velocity[i])
        
        tidenames = ['Neap','Intermediate','Spring']
        for i in range(3):
            if thetisfilename == 'vsotf_otfmesh':
                axs[i].plot( otf_time[i], thetis_velocity3[51*i:51*(i+1)],  label='VSotf')
            else:
                axs[i].plot( data_time[i], thetis_velocity3[151*i:151*(i+1)],  label='Thetis')
            axs[2].set_xlabel('$Time$ $[date]$', fontsize=18)
            axs[i].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
            axs[i].set_xticks(data_time[i][0:-1:48])
            yy=np.linspace(0,2.0,5)
            axs[i].set_yticks(yy)
            #axs[i].set_ylim((0.00,2.00))
            axs[i].set_ylabel('$u$ $[m/s]$', fontsize=18)
            axs[i].xaxis.set_tick_params(labelsize=15)
            axs[i].yaxis.set_tick_params(labelsize=15)
            axs[i].set_title('Velocity magnitude time series at B2 -- '+tidenames[i]+' tide',fontsize=18)
            axs[0].legend(loc='best',fontsize=14.5,ncol=4)
        plt.savefig('sediment-velocity.png', dpi =300)
#axs.plot( data_time[3], thetis_velocity,  label='Thetis')

show()


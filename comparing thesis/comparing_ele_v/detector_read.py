from cProfile import label
import datetime
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import h5py
import numpy as np
import pandas as pd
import matplotlib.dates as mdates

zw = matplotlib.font_manager.FontProperties(fname='/media/can/can_disk/simsun.ttc',size=15)
# yw = matplotlib.font_manager.FontProperties(fname='/media/can/can_disk/times.ttf')
plt.rc('font',family='Times New Roman')

time60min=[]

startneap = datetime.datetime(2013,8,16,10,0,0)#start time
stopneap = datetime.datetime(2013,8,17,11,0,0)#end time
delta = datetime.timedelta(minutes=60)
dateneap = mpl.dates.drange(startneap,stopneap,delta)
time60min.append(dateneap)

startmeso = datetime.datetime(2013,8,19,14,0,0)
stopmeso = datetime.datetime(2013,8,20,15,0,1)
delta = datetime.timedelta(minutes=60)
datemeso = mpl.dates.drange(startmeso,stopmeso,delta)
time60min.append(datemeso)

startspring = datetime.datetime(2013,8,23,10,0,0)
stopspring = datetime.datetime(2013,8,24,11,0,1)
delta = datetime.timedelta(minutes=60)
dataspring = mpl.dates.drange(startspring,stopspring,delta)
time60min.append(dataspring)

time5min = []

startneap = datetime.datetime(2013,8,16,10,0,0)#start time
stopneap = datetime.datetime(2013,8,17,11,0,0)#end time
delta = datetime.timedelta(minutes=5)
dateneap = mpl.dates.drange(startneap,stopneap,delta)
time5min.append(dateneap)

startmeso = datetime.datetime(2013,8,19,14,0,0)
stopmeso = datetime.datetime(2013,8,20,15,0,1)
delta = datetime.timedelta(minutes=5)
datemeso = mpl.dates.drange(startmeso,stopmeso,delta)
time5min.append(datemeso)

startspring = datetime.datetime(2013,8,23,10,0,0)
stopspring = datetime.datetime(2013,8,24,11,0,1)
delta = datetime.timedelta(minutes=5)
dataspring = mpl.dates.drange(startspring,stopspring,delta)
time5min.append(dataspring)

A = ['x','z','d']
B = ['n','i','s']
C = ['小','中','大']
D = ['0.5','1.0','2.0']

domainnames = ['A','B','C']


stationnames=['B1','B3']#,'B3','B4','B5']#,'A1','A2','A3','A4','A5']
labelnames = ['C1','C2']
for ln,sn in enumerate(stationnames):
    fig,axs = plt.subplots(3,1,figsize=(16,9))
    axs=axs.reshape(-1)
    fig.tight_layout(pad=6,h_pad=3.5)
    
    for aa,bb,cc,dd in zip(A,B,C,D):
        ii = A.index(aa)
        ax = axs[ii]
        # ax1 = ax.twinx()
        df_file = pd.read_excel('./sediment.xls',header=7,usecols='J',sheet_name=sn+aa)
        sed_m = np.array(df_file['垂向平均'][1:])
        ax.scatter(time60min[ii],sed_m,s=20,c='k',label='Measurement')
        c_m = mean(sed_m)
        for domainname in domainnames:
            det_file = h5py.File('../../outputs/0.validation/sediment-nis/'+domainname+'-'+bb+'/diagnostic_detectors.hdf5')
            sediment,vx,vy = [],[],[]
            for name, data in det_file.items():
                # print(name)
                if name == sn:
                    sediment.append(data[:574*6,3])
            # det_file = '../../../0.validation/continuous-sediment-exner2/diagnostic_detectors.hdf5'
            # df = h5py.File(det_file,'r+')
            # for name, data in df.items():
            #     if name == sn:
            #         sediment[0]=np.append(sediment[0],data[:,3])

            sed1 = list(sediment[0][6*6:6*6+301])
            ax.plot(time5min[ii],sed1,label='计算域'+domainname)
            c_s = mean(sed1)


            axs[2].set_xlabel('日期', fontsize=18,fontproperties=zw)
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
            ax.set_xticks(time60min[ii][0:-1:6])
            if ii == 2:
                yy= np.arange(0.0,2.51,0.50)
            else:
                yy= np.arange(0.0,2.51,0.25)
            ax.set_yticks(yy)
            ax.set_ylim(-0.05,float(dd))
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
            ax.set_ylabel('含沙量 $[kg/m^3]$', fontsize=18,fontproperties=zw)
            ax.xaxis.set_tick_params(labelsize=15)
            ax.yaxis.set_tick_params(labelsize=15)
            ax.set_title('含沙量测点'+labelnames[ln]+'处'+cc+'潮期间含沙量变化',fontsize=18,fontproperties=zw)
            axs[0].legend(loc='best',prop=zw,ncol=4)
            print(domainname+'\t',cc+'潮'+labelnames[ln]+'处：{0:.2f}'.format(abs(c_m-c_s)/c_m))
    plt.savefig(labelnames[ln]+'.png',dpi=300)

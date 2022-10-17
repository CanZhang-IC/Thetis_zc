import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score

zw = matplotlib.font_manager.FontProperties(fname='/media/can/can_disk/simsun.ttc',size=25)
# yw = matplotlib.font_manager.FontProperties(fname='/media/can/can_disk/times.ttf')
plt.rc('font',family='Times New Roman')

axis_size, axis_label_size, title_size, legend_size =  20, 20, 20, 18

for angle1 in range(0,31,10):
    fig,axs = plt.subplots(1,6,figsize=(16,9))
    axs = axs.reshape(-1)
    fig.tight_layout(pad = 5.5, h_pad = 3)
    fig.subplots_adjust(wspace=0.0,hspace=0.02)
    xx_sticks=np.linspace(0,1.2,5)
    yy_sticks=np.linspace(-1,1,5)
    yy_number = 15

    locations=[]
    names=[]
    xx=[1,2,4,6,8,10]
    yy=[-1,-0.8,-0.6,-0.5,-0.4,-0.2,-0.1, 0, 0.1,0.2,0.4,0.5,0.6,0.8,1]

    for i in xx:
        for j in yy:
            locations.append((25+0.3*i,2.5+0.3*j))
            names.append('name'+str((i,j)))

    ###read measured data
    measure_velocity = []
    #combine data from 1,2,3,4,6,8,10D
    for x_distance in xx:
        # print(x_distance)
        measure_read = pd.read_excel('./measurement/one turbine/'+str(angle1)+'/'+str(x_distance)+'D.xls')
        measure_array = np.array(measure_read)
        for i in measure_array:
            measure_velocity.append(i[1]/0.33)

    for i in range(6):
        axs[i].plot(measure_velocity[yy_number*i:yy_number*(i+1)],yy, 'ko',  label = '测量结果')


    ###Read Thetis data###
    def_dirfile = h5py.File('../../../outputs/0.validation/experiment_turbine/0.6-1.5structure/'+str(float(angle1))+'/diagnostic_detectors.hdf5','r+')

    xvelocity = []
    yvelocity = []      
    for name, data in def_dirfile.items():
        if name == 'time':
            pass
        else:
            xvelocity.append(data[-1][1])
            yvelocity.append(data[-1][2])     
    thetis_velocity_all=[]
    for i in range(len(xvelocity)):
        thetis_velocity_all.append(xvelocity[i]/0.33)      

    thetis_velocity=[]
    #devide data in 1,2,4,6,8,10D
    for i in range(6):
        #velocity in original order
        original_x = thetis_velocity_all[yy_number*i:yy_number*(i+1)]
        #velocity in reversed order
        reversed_x=[]
        for ii in reversed(original_x):
            reversed_x.append(ii)
        # combine half original order and half reversed order.
        for iii in range(int(yy_number/2)+1):
            thetis_velocity.append(reversed_x[iii])
        for iiii in range(int(yy_number/2)):
            thetis_velocity.append(original_x[iiii])

    v_10d = thetis_velocity[15:30]
    del thetis_velocity[15:30]
    v_12d = thetis_velocity[15:30]
    del thetis_velocity[15:30]
    for i in v_10d:
        thetis_velocity.append(i)
    for i in v_12d:
        thetis_velocity.append(i)


    for i in range(6):
        axs[i].plot(thetis_velocity[yy_number*i:yy_number*(i+1)], yy ,'k', label = '数值结果')


        #设置坐标轴坐标、范围以及字体大小
        axs[i].set_xticks([0,0.5,1.0])
        axs[i].set_yticks([-1.0,-0.5,0,0.5,1.0])
        axs[i].set_xlim(0,1.5)
        axs[i].set_ylim(-1.05,1.05)
        axs[i].xaxis.set_tick_params(labelsize=axis_size)
        axs[i].yaxis.set_tick_params(labelsize=axis_size)
        
        #设置标题以及字体大小
        axs[i].set_title(str(xx[i])+'D',fontsize=title_size,x=0.18,y=0.93)

        if i == 0:
            axs[i].set_ylabel('$Y/D$',fontsize=axis_label_size)
            axs[i].set_yticks([-1.0,-0.5,0,0.5,1.0])
            axs[i].set_xticks([0,0.5,1.0])
        else:
            axs[i].yaxis.set_tick_params(labelleft=False) #不显示y轴上的刻度值，但是可以保留网格
        ###计算RMSE和R2###
        thetis_rmse_elevation1=mean_squared_error(np.array(thetis_velocity[yy_number*i:yy_number*(i+1)]),np.array(measure_velocity[yy_number*i:yy_number*(i+1)]), squared=False)
        thetis_r2_elevation1=r2_score(np.array(thetis_velocity[yy_number*i:yy_number*(i+1)]),np.array(measure_velocity[yy_number*i:yy_number*(i+1)]))
        print(angle1,xx[i],thetis_rmse_elevation1,thetis_r2_elevation1)

    axs[2].set_xlabel('$u/u_0$',fontsize=axis_label_size,x = 1) #设置坐标轴名称及位置
    axs[5].set_xticks([0,0.5,1.0,1.5])
    plt.legend( bbox_to_anchor=(-3.2, 1),#设置图例的位置
    loc=3, borderaxespad=0,
        ncol=4,#水平放置图例的数量
        frameon = False, #不显示图例边框
        prop=zw
        )
    plt.savefig('./'+str(angle1)+'.jpg',dpi=300)
    



    
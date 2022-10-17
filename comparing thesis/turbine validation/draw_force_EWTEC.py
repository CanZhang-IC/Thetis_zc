import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import h5py
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score

zw = matplotlib.font_manager.FontProperties(fname='/media/can/can_disk/simsun.ttc',size=15)
# yw = matplotlib.font_manager.FontProperties(fname='/media/can/can_disk/times.ttf')
plt.rc('font',family='Times New Roman')

fig,axs = plt.subplots(1,1,figsize=(9,9))
fig.subplots_adjust(wspace=0.3)

angle = [0,10,20,30]

### analytical
D,Dt,u0 = 0.3,0.6,0.33
smooth_angles=range(31)
fx_a_plot = [0.5*np.pi*(D/2)**2*Dt*(u0*np.cos(a/180*np.pi))**2*1000 for a in smooth_angles]
fy_a_plot = [0.5*np.pi*(D/2)**2*Dt*(u0*np.sin(a/180*np.pi))**2*1000 for a in smooth_angles]
axs[0].plot(smooth_angles,fx_a_plot,'k-',label='解析结果：x方向分力')
# axs[0].plot(smooth_angles,fy_a_plot,'k-',label='解析结果：y方向分力')

fx_a = [0.5*np.pi*(D/2)**2*Dt*(u0*np.cos(a/180*np.pi))**2*1000 for a in angle]
fy_a = [0.5*np.pi*(D/2)**2*Dt*(u0*np.sin(a/180*np.pi))**2*1000 for a in angle]


### thetis result
fx_t = [2.228401286,2.265867075,2.128419305,1.934036784] 
fy_t = [0,0.3995335,0.774681273,1.116616658]
axs[0].scatter(angle,fx_t,s = 60,c='k',marker='s',label='数值结果：x方向分力')
# axs[0].scatter(angle,fy_t,s = 60,c='k',marker='s',label='数值结果：y方向分力')

### measurement
fx_m = [2.255957004,2.25052013,2.096790034,1.672413991] 
fy_m = [0, 0.70499143, 1.399274461, 1.974595044]
axs[0].scatter(angle,fx_m,s = 60,c='k',marker='P',label='测量结果：x方向分力')
# axs[0].scatter(angle,fy_t,s = 60,c='k',marker='P',label='测量结果：y方向分力')

axs[0].set_ylabel('水轮机的受力 $[N]$',fontsize=20,fontproperties=zw)
axs[0].set_xlabel('偏航角度 $[\circ]$',fontsize=20,fontproperties=zw)



# #### EWTEC图 ###

# ### analytical
# D,Dt,u0 = 20,0.6,2
# smooth_angles=range(0,91,5)
# fx_a_plot = [0.5*np.pi*(D/2)**2*Dt*(u0*np.cos(a/180*np.pi))**2 for a in smooth_angles]
# # print(fx_a_plot)
# axs[1].plot(smooth_angles,fx_a_plot,'k-',label='解析结果：x方向分力')

# power_original = []
# for angle in range(0,81,10):
#     file_dir = '../../../outputs/4.yaw/Yaw_Ideal/Conference/conference_mesh2_5-10/f30-cos22/angle'+str(angle)+'H'+str(40)
#     file_value = h5py.File(file_dir+'/diagnostic_turbine.hdf5')
#     for name,value in file_value.items():
#         # print(name)
#         if name == 'F_force':
#             power_original.append(-float(value[-1]))
# theangle = range(0,91,10)
# axs[1].scatter(theangle,power_original+[0],s = 60,c='k',marker='s',label='数值结果：x方向分力')

# axs[1].set_ylabel('水轮机的受力 $[kN]$',fontsize=20,fontproperties=zw)
# axs[1].set_xlabel('偏航角度 $[\circ]$',fontsize=20,fontproperties=zw)
#设置坐标轴坐标、范围以及字体大小
for i in range(1):
    # axs[i].set_xticks([0,0.5,1.0])
    # axs[i].set_yticks([-1.0,-0.5,0,0.5,1.0])
    # axs[i].set_xlim(0,1.5)
    # axs[i].set_ylim(-1.05,1.05)
    axs[i].xaxis.set_tick_params(labelsize=20)
    axs[i].yaxis.set_tick_params(labelsize=20)
    axs[i].legend(prop=zw)

plt.savefig('./受力.png',dpi=300)

###计算RMSE和R2###
thetis_rmse_elevation1=mean_squared_error(np.array(fx_t),np.array(fx_m), squared=False)
thetis_rmse_elevation2=mean_squared_error(np.array(fx_t),np.array(fx_a), squared=False)

# thetis_r2_elevation1=r2_score(np.array(thetis_velocity[yy_number*i:yy_number*(i+1)]),np.array(measure_velocity[yy_number*i:yy_number*(i+1)]))
print('数值vs测量',thetis_rmse_elevation1)
print('数值vs解析',thetis_rmse_elevation2)





# D,Dt,u0 = 20,0.6,2
# fx_a_EWTEC = [0.5*np.pi*(D/2)**2*Dt*(u0*np.cos(a/180*np.pi))**2 for a in theangle]
# thetis_rmse_elevation3=mean_squared_error(np.array(power_original+[0]),np.array(fx_a_EWTEC), squared=False)
# print('EWTEC：数值vs解析',thetis_rmse_elevation3)

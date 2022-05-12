
import datetime
from pylab import *
import matplotlib
import h5py
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from math import atan
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
import pandas as pd


#read measured data
file_name ='B5.csv'
velocity  = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(1),skiprows=1) 
direction = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(2),skiprows=1) 
# data formatting:
measured_velocity = []
measured_direction = []
# remove empty values
for i in range(len(velocity)):
    if velocity[i] != '-':
        measured_velocity.append(float(velocity[i]))
        measured_direction.append(float(direction[i]))


# ###read OTF data
# file_name = 'v_otf.csv'
# otf_v = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(2),skiprows=0) 
# otf_d = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(3),skiprows=0) 
# # data formatting:
# all_otf_v = []
# all_otf_d = []
# # remove empty values
# for i in range(len(otf_v)):
#     if otf_v[i] != '-':
#         all_otf_v.append(float(otf_v[i]))
#         all_otf_d.append(float(otf_d[i]))
# otf_velocity_magnitude=[]
# otf_velocity_direction=[]
# for i in range(149,200):
#     otf_velocity_magnitude.append(all_otf_v[i])
#     otf_velocity_direction.append(all_otf_d[i])
# for i in range(300,351):
#     otf_velocity_magnitude.append(all_otf_v[i])
#     otf_velocity_direction.append(all_otf_d[i])
# for i in range(484,535):
#     otf_velocity_magnitude.append(all_otf_v[i])
#     otf_velocity_direction.append(all_otf_d[i])

# ### Far-field B5 Read OTF data ###
# file_name = 'v_otf.csv'
# otf_xv = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(5),skiprows=0) 
# otf_yv = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(6),skiprows=0)
# # data formatting:
# all_otf_v = []
# all_otf_d = []
# # remove empty values
# for i in range(len(otf_xv)):
#     all_otf_v.append(sqrt(float(otf_xv[i])**2+float(otf_yv[i])**2))
#     one_direction = np.arctan2(float(otf_xv[i]),float(otf_yv[i]))* 180 / np.pi
#     if one_direction < 0:
#         all_otf_d.append(one_direction+360)
#     else:
#         all_otf_d.append(one_direction)

# otf_velocity_magnitude=[]
# otf_velocity_direction=[]
# for i in range(147,198):
#     otf_velocity_magnitude.append(all_otf_v[i])
#     otf_velocity_direction.append(all_otf_d[i])
# for i in range(298,349):
#     otf_velocity_magnitude.append(all_otf_v[i])
#     otf_velocity_direction.append(all_otf_d[i])
# for i in range(482,533):
#     otf_velocity_magnitude.append(all_otf_v[i])
#     otf_velocity_direction.append(all_otf_d[i])

### Far-field A4 Read OTF data ###
file_name = 'v_otf.csv'
otf_xv = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(8),skiprows=0) 
otf_yv = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(9),skiprows=0)
# data formatting:
all_otf_v = []
all_otf_d = []
# remove empty values
for i in range(len(otf_xv)):
    all_otf_v.append(sqrt(float(otf_xv[i])**2+float(otf_yv[i])**2))
    one_direction = np.arctan2(float(otf_xv[i]),float(otf_yv[i]))* 180 / np.pi
    if one_direction < 0:
        all_otf_d.append(one_direction+360)
    else:
        all_otf_d.append(one_direction)

otf_velocity_magnitude=[]
otf_velocity_direction=[]
for i in range(147,198):
    otf_velocity_magnitude.append(all_otf_v[i])
    otf_velocity_direction.append(all_otf_d[i])
for i in range(298,349):
    otf_velocity_magnitude.append(all_otf_v[i])
    otf_velocity_direction.append(all_otf_d[i])
for i in range(482,533):
    otf_velocity_magnitude.append(all_otf_v[i])
    otf_velocity_direction.append(all_otf_d[i])

otf_rmse_velocity=[]
otf_rmse_direction=[]
#neap,meso,spring
fig, axs = plt.subplots(3, 1, figsize=(16,9))
axs = axs.reshape(-1) # reshape so that we can iterate below over axs[i] instead of ax[i,j]
fig.tight_layout(pad=3,h_pad=1.5)
for i in range(3):
    otf_rmse_velocity.append(np.sqrt(mean_squared_error(np.array(measured_velocity[151*i:151*(i+1):3]),np.array(otf_velocity_magnitude[51*i:51*(i+1)]))))
    #otf_rmse_velocity.append(r2_score(np.array(measured_velocity[151*i:151*(i+1):3]),np.array(otf_velocity_magnitude[51*i:51*(i+1)])))
    count1 = [i+1 for i in range(51)]
    axs[i].plot(count1,measured_velocity[151*i:151*(i+1):3],'ko')
    axs[i].plot(count1,otf_velocity_magnitude[51*i:51*(i+1)])
    direction_diff = (np.array(measured_direction[151*i:151*(i+1):3]) - np.array(otf_velocity_direction[51*i:51*(i+1)])+180) % 360 -180
    otf_rmse_direction.append(np.sqrt(np.sum(direction_diff**2)/direction_diff.shape[0]))
#all
otf_rmse_velocity.append(0)#np.sqrt(mean_squared_error(np.array(velocity_magnitude[::3]),np.array(otf_velocity_magnitude))))
otf_rmse_direction.append(0)#np.sqrt(mean_squared_error(np.array(measured_direction)/360,np.array(otf_velocity_direction)/360)))
print('The OpenTidalFarm simulated rmse of velocity  during neap tide is {0:.3f}, during meso tide is {1:.3f}, during spring tide is {2:.3f} and totally is {3:.3f}.'.
    format(otf_rmse_velocity[0],otf_rmse_velocity[1],otf_rmse_velocity[2],otf_rmse_velocity[3]))
print('The OpenTidalFarm simulated rmse of direction during neap tide is {0:.3f}, during meso tide is {1:.3f}, during spring tide is {2:.3f} and totally is {3:.3f}.'.
    format(otf_rmse_direction[0],otf_rmse_direction[1],otf_rmse_direction[2],otf_rmse_direction[3]))

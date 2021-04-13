import h5py
import matplotlib.pyplot as plt
import numpy as np

# a,b = 16,9
# fig, axs = plt.subplots(1, 1, figsize=(a,b))

# file_dir = '../../../outputs/4.yaw/Yaw_Ideal/op-conference_mesh2-5_40/f30-cos00/optimisation-aligned-angle'
# file_value = h5py.File(file_dir+'/diagnostic_controls.hdf5')

# for name,value in file_value.items():
#     angles = value[-1]
# print(list(angles))
# # angles = [0] *12
# turbine_location = []
# for i in range(850,1200,100):
#     for j in range(250, 400, 50):
#         turbine_location.append([i,j])
# # for i in range(850,1200,200):
# #     for j in range(250, 400, 50):
# #         turbine_location.append([i,j])
# # for i in range(950,1200,200):
# #     for j in range(275, 400, 50):
# #         turbine_location.append([i,j])

# for i in range(len(angles)):
#     point1_x = turbine_location[i][0] + 15*np.cos((90+angles[i])/180*np.pi)
#     point2_x = turbine_location[i][0] + 15*np.cos((90+angles[i]+180)/180*np.pi)
#     point1_y = turbine_location[i][1] + 15*np.sin((90+angles[i])/180*np.pi)
#     point2_y = turbine_location[i][1] + 15*np.sin((90+angles[i]+180)/180*np.pi)
#     axs.plot(turbine_location[i][0],turbine_location[i][1],'o')
#     axs.plot([point1_x,point2_x],[point1_y,point2_y],'k',linewidth =12)

# axs.legend(loc='best',fontsize=14.5)
# axs.set_xlabel('X ($m$)', fontsize=18)
# axs.xaxis.set_tick_params(labelsize=15)
# axs.set_ylabel('Y ($m$)', fontsize=18)
# axs.yaxis.set_tick_params(labelsize=15)
# # xx=np.linspace(0,2000,100)
# # axs.set_xticks(xx)
# # yy=np.linspace(0,600,50)
# # axs.set_xticks(yy)
# plt.savefig('./cos00/aligned-angle.jpg',dpi=300)

a,b = 16,9
fig, axs = plt.subplots(1, 1, figsize=(a,b))

file_dir = '../../../outputs/4.yaw/Yaw_Ideal/op-conference_mesh2-5_40/f30-cos00/optimisation-staggered-both'
file_value = h5py.File(file_dir+'/diagnostic_controls.hdf5')

for name,value in file_value.items():
    location_angle = value[-1]
locations = list(location_angle[:24])

angles = list(location_angle[24:]) 

print(angles,locations)
turbine_location = []
for i in range(int(len(locations)/2)):
    turbine_location.append([locations[2*i],locations[2*i+1]])


for i in range(len(angles)):
    point1_x = turbine_location[i][0] + 15*np.cos((90+angles[i])/180*np.pi)
    point2_x = turbine_location[i][0] + 15*np.cos((90+angles[i]+180)/180*np.pi)
    point1_y = turbine_location[i][1] + 15*np.sin((90+angles[i])/180*np.pi)
    point2_y = turbine_location[i][1] + 15*np.sin((90+angles[i]+180)/180*np.pi)
    axs.plot(turbine_location[i][0],turbine_location[i][1],'o')
    axs.plot([point1_x,point2_x],[point1_y,point2_y],'k',linewidth =12)


axs.set_xlabel('X ($m$)', fontsize=18)
axs.xaxis.set_tick_params(labelsize=15)
axs.set_ylabel('Y ($m$)', fontsize=18)
axs.yaxis.set_tick_params(labelsize=15)
# xx=np.linspace(0,2000,100)
# axs.set_xticks(xx)
# yy=np.linspace(0,600,50)
# axs.set_xticks(yy)
plt.savefig('./cos00/staggered_both_lines.jpg',dpi=300)
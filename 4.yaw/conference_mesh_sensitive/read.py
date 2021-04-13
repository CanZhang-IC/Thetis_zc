import h5py
import matplotlib.pyplot as plt
import numpy as np

mesh_list1 = ['conference_mesh2_2-10','conference_mesh2_5-10','conference_mesh2_5-20','conference_mesh2_5-40','conference_mesh2_10-40']

# fig, axs = plt.subplots(1, 1, figsize=(12,9))
# # axs = axs.reshape(-1) # reshape so that we can iterate below over axs[i] instead of ax[i,j]
# # fig.tight_layout(pad=4,h_pad=3.5)
# axs.plot(range(0,91,5),[615.4212368,
# 609.7614534,
# 592.889474,
# 565.1678046,
# 527.3139539,
# 480.5085821,
# 426.4581115,
# 367.3770433,
# 305.8809087,
# 244.8019649,
# 186.9529437,
# 134.8703598,
# 90.56996951,
# 55.34433908,
# 29.62695229,
# 12.9395357,
# 3.930013159,
# 0.498601348,
# 1.73103E-46
# ],'o',label='analytical')

# print('\t','\t','Original')
# for ii,mesh_r in enumerate(mesh_list1):
#     mesh_label = str(mesh_r)
#     print('mesh:',mesh_label)
#     for H in range(40,41,10):
#         power_original = []
#         for angle in range(0,81,10):
#             file_dir = '../../../outputs/4.yaw/Yaw_Ideal/'+str(mesh_r)+'/f30-cos22/angle'+str(angle)+'H'+str(H)
#             file_value = h5py.File(file_dir+'/diagnostic_turbine.hdf5')
#             for name,value in file_value.items():
#                 # print(name)
#                 if name == 'current_power':
#                     power_original.append(float(value[-1]))
#                     print('Angle:', angle, '\t', 'H:', H, '\t',float(value[-1]))
#         theangle = range(0,91,10)
#         axs.plot(theangle,power_original+[0],label='Mesh'+str(ii+1))
        


# axs.legend(loc='best',fontsize=14.5)
# axs.set_xlabel('Yaw angle ($\circ$)', fontsize=18)
# axs.xaxis.set_tick_params(labelsize=15)
# axs.set_ylabel('Power output ($kW$)', fontsize=18)
# axs.yaxis.set_tick_params(labelsize=15)
# xx=np.linspace(0,90,10)
# axs.set_xticks(xx)
# # axs.set_title('Numerical power output under different yaw condition', fontsize=18)

# plt.savefig('power_mesh.jpg', dpi = 300)
# # plt.show()

fig, axs = plt.subplots(1, 1, figsize=(12,9))
# axs = axs.reshape(-1) # reshape so that we can iterate below over axs[i] instead of ax[i,j]
# fig.tight_layout(pad=4,h_pad=3.5)
analytical_force = [376.9911184,
374.1274473,
365.6234453,
351.737502,
332.8915349,
309.6581692,
282.7433388,
252.9648374,
221.2274696,
188.4955592,
155.7636489,
124.026281,
94.24777961,
67.33294927,
44.09958353,
25.25361643,
11.36767317,
2.863671092,
1.41349E-30
]
axs.plot(range(0,91,5), analytical_force,'o',label='Analytical')

print('\t','\t','Original')
for ii,mesh_r in enumerate(mesh_list1):
    mesh_label = str(mesh_r)
    print('mesh:',mesh_label)
    for H in range(40,41,10):
        force_original = []
        for angle in range(0,81,10):
            file_dir = '../../../outputs/4.yaw/Yaw_Ideal/'+str(mesh_r)+'/f30-cos22/angle'+str(angle)+'H'+str(H)
            file_value = h5py.File(file_dir+'/diagnostic_turbine.hdf5')
            for name,value in file_value.items():
                # print(name)
                if name == 'F_force':
                    force_original.append(-float(value[-1]))
                    print('Angle:', angle, '\t', 'H:', H, '\t',float(value[-1]))
        theangle = range(0,91,10)
        force_original = force_original + [0]
        axs.plot(theangle,force_original,label='Mesh'+str(ii+1))
        Rmse_results = np.sqrt(sum([(abs(a)-abs(b))**2 for a,b in zip(analytical_force[::2],force_original)])/len(theangle))
        print(Rmse_results)
        


axs.legend(loc='best',fontsize=14.5)
axs.set_xlabel('Yaw angle ($\circ$)', fontsize=18)
axs.xaxis.set_tick_params(labelsize=15)
axs.set_ylabel('Turbine thrust ($N$)', fontsize=18)
axs.yaxis.set_tick_params(labelsize=15)
xx=np.linspace(0,90,10)
axs.set_xticks(xx)
# axs.set_title('Numerical turbine thrust under different yaw condition', fontsize=18)

# plt.savefig('force_mesh.jpg', dpi = 300)
# plt.show()
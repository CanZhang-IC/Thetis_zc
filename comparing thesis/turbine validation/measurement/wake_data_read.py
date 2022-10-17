import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error

fig,axs = plt.subplots(2,4,figsize=(10,10))
axs = axs.reshape(-1)
fig.tight_layout(pad = 3, h_pad = 1.5)
xx_sticks=np.linspace(0,1.2,5)
yy_sticks=np.linspace(-1,1,5)
yy_number = 15

for angle in ['0 0']:#, '0 10', '0 20', '0 30', '10 0', '20 0', '30 0']:
    print(angle)

    xx=[1,2,3,4,6,8,10]
    yy=[-1,-0.8,-0.6,-0.5,-0.4,-0.2,-0.1, 0, 0.1,0.2,0.4,0.5,0.6,0.8,1]

    ###read measured data
    #For two turbines. Combine data from 1,2,3,4,6,8,10D
    measure_velocity = []
    for x_distance in xx:
        measure_read = pd.read_excel('./two turbines/'+str(angle)+'/'+str(x_distance)+'D.xls')
        measure_array = np.array(measure_read)
        for i in measure_array:
            measure_velocity.append(i[1]/0.33)

    for i in range(7):
        axs[i].plot(measure_velocity[yy_number*i:yy_number*(i+1)],yy,  label = 'Measure_'+str(angle))
        axs[i].set_xlabel(str(xx[i])+'D', fontsize = 14)
        if i == 0 or i ==5:
            axs[i].set_ylabel('$L/L_0$',fontsize = 14)
        axs[i].set_xticks(xx_sticks)
        axs[i].set_yticks(yy_sticks)
        axs[i].legend(loc = 'best', fontsize = 6)

    #For one turbine. Combine data from 1,2,3,4,6,8,10D
    measure_velocity = []
    for x_distance in xx:
        # print(x_distance)
        measure_read = pd.read_excel('./one turbine/'+str(angle[0])+'/'+str(x_distance)+'D.xls')
        measure_array = np.array(measure_read)
        for i in measure_array:
            measure_velocity.append(i[1]/0.33)

    for i in range(7):
        axs[i].plot(measure_velocity[yy_number*i:yy_number*(i+1)],yy,  label = 'Measure_'+str(angle[0]))
        axs[i].set_xlabel(str(xx[i])+'D', fontsize = 14)
        if i == 0 or i ==5:
            axs[i].set_ylabel('$L/L_0$',fontsize = 14)
        axs[i].set_xticks(xx_sticks)
        axs[i].set_yticks(yy_sticks)
        axs[i].legend(loc = 'best', fontsize = 6)
  
    plt.savefig('measure.jpg',dpi=300)



    
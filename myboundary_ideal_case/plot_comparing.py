import datetime
from pylab import *
import matplotlib
import h5py
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from math import sqrt,atan



      
det_file = "outputs/diagnostic_detectors.hdf5"
df = h5py.File(det_file, 'r+')
xvelocity=[]
yvelocity=[]
for name, data in df.items():
    if name == 'point':
        xvelocity.append(data[:,1])
        yvelocity.append(data[:,2])
smalldomain_velocity=[]

for i in range(len(xvelocity[0])):
    smalldomain_velocity.append(sqrt(xvelocity[0][i]**2+yvelocity[0][i]**2))    

det_file = "large_domain/outputs/diagnostic_detectors.hdf5"
df = h5py.File(det_file, 'r+')
xvelocity=[]
yvelocity=[]
for name, data in df.items():
    if name == 'point':
        xvelocity.append(data[:,1])
        yvelocity.append(data[:,2])
largedomain_velocity=[]

for i in range(len(xvelocity[0])):
    largedomain_velocity.append(sqrt(xvelocity[0][i]**2+yvelocity[0][i]**2)) 

t = [i for i in range(len(smalldomain_velocity))]
plt.plot(t, smalldomain_velocity, 'b', label='small_domain-result')
plt.plot(t, largedomain_velocity, 'r.', label='large_domain-result')
plt.legend(loc='best')
show()


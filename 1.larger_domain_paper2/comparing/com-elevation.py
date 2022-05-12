from uptide import *
import datetime
import math
from pylab import *
import pandas as pd

import matplotlib
import h5py
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
import shutil

#read measured data
file_name = 'B.csv'
water_level = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(2),skiprows=1) # water level (cm) in 2nd column
time = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(1),skiprows=1) # time in seconds from initial time in 1rd column
# data formatting:
both = []
# remove empty values
for i in range(len(time)):
    if water_level[i] != '-':
        both.append([water_level[i], time[i]])       
x = []
t = []
for i in range(len(both)):
    x.append(float(both[i][0])) # converting water level from cm to m
    t.append(float(both[i][1]))
plot(t,x,'ko')

#read Thetis data
thetisfilename='onemin-8cores
  
det_file = "../outputs/"+thetisfilename+"/diagnostic_detectors.hdf5"
df = h5py.File(det_file, 'r+')

dt=60
spin = int(3*24*60*60/dt) # 3 days as simulation begin at 09/08/2013 00:00, measure data begin at 12/08/2013 00:00
t_simu = np.arange(28800,1324800,60)
# 8 hours for different time zone [0,1296000,60] to [28800,1324800,60]. 1296000 is the final time of the measured data standing for 27/08/2013 00:00
waterlevel_thetis=[]
for name, data in df.items():
	if name == 'B':
		waterlevel_thetis.append(data[spin:,0])
plot(t_simu,waterlevel_thetis[0],'b')

#read OTF data
file_name = 'otf_water_level.csv'
otf_water_level = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(1),skiprows=1) # water level (cm) in 2nd column
otf_time = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(0),skiprows=1) # time in seconds from initial time in 1rd column
# data formatting:
otf_both = []
# remove empty values
for i in range(len(otf_time)):
    if otf_water_level[i] != '-':
        otf_both.append([otf_water_level[i], otf_time[i]])       
otf_x = []
otf_t = []
for i in range(len(otf_both)):
    otf_x.append(float(otf_both[i][0])) # converting water level from cm to m
    otf_t.append(float(otf_both[i][1]))
plot(otf_t,otf_x,'g')

show()


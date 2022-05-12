
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

rmse_velocity=[]
rmse_direction=[]

names=[
	'onemin-8cores-huluthreepart'
      ]
plotnames=['test1']
for filename in names:
	stationnames=['B5']
	for stationname in stationnames:
		print(stationname)
		#read measured data
		file_name =stationname+'.csv'
		velocity = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(1),skiprows=1) # velocity in 2nd column
		direction = np.loadtxt(file_name, dtype=str, delimiter=',', usecols=(2),skiprows=1)

		# data formatting:
		both = []
		# remove empty values
		for i in range(len(velocity)):
		    if velocity[i] != '-':
		        both.append([velocity[i], direction[i]])

		velocity_magnitude = []
		velocity_direction= []

		for i in range(len(both)):
		    velocity_magnitude.append(float(both[i][0]))
		    velocity_direction.append(float(both[i][1]))

		det_file = "../outputs/"+filename+"/diagnostic_detectors.hdf5"
		df = h5py.File(det_file, 'r+')
		xvelocity=[]
		yvelocity=[]
		for name, data in df.items():
		    if name == stationname:       
		        xvelocity.append(data[:,1])
		        yvelocity.append(data[:,2])

		thetis_velocity=[]
		thetis_direction=[]
		for i in range(len(xvelocity[0])):
		    thetis_velocity.append(sqrt(xvelocity[0][i]**2+yvelocity[0][i]**2))
		    if yvelocity[0][i] < 0 :
		        thetis_direction.append(atan(xvelocity[0][i]/yvelocity[0][i])/np.pi*180+180)
		    else:
		        if xvelocity[0][i] < 0:
		            thetis_direction.append(atan(xvelocity[0][i]/yvelocity[0][i])/np.pi*180+360)
		        else:
		            thetis_direction.append(atan(xvelocity[0][i]/yvelocity[0][i])/np.pi*180)
		
		thetis_velocity3=[]
		thetis_direction3=[]
		if filename=='minimumdepth(2)-manning(0.02-0.2)-viscosity(5-1000)-simulation-test3-1min':
			for i in range(10200,11701,10):
				thetis_velocity3.append(np.mean(thetis_velocity[i:i+10]))
				thetis_direction3.append(np.mean(thetis_direction[i:i+10]))
			for i in range(14760,16261,10):
				thetis_velocity3.append(np.mean(thetis_velocity[i:i+10]))
				thetis_direction3.append(np.mean(thetis_direction[i:i+10]))
			for i in range(20280,21781,10):
				thetis_velocity3.append(np.mean(thetis_velocity[i:i+10]))
				thetis_direction3.append(np.mean(thetis_direction[i:i+10]))
		else:
			for i in range(2040,2341,2):
				thetis_velocity3.append(np.mean(thetis_velocity[i:i+2]))
				thetis_direction3.append(np.mean(thetis_direction[i:i+2]))
			for i in range(2952,3253,2):
				thetis_velocity3.append(np.mean(thetis_velocity[i:i+2]))
				thetis_direction3.append(np.mean(thetis_direction[i:i+2]))
			for i in range(4056,4357,2):
				thetis_velocity3.append(np.mean(thetis_velocity[i:i+2]))
				thetis_direction3.append(np.mean(thetis_direction[i:i+2]))

		'''
		for i in range(len(thetis_direction3)):
			if i == 26 :
				thetis_direction3[i] = 150 + thetis_direction3[i]
			if i == 28 :
				thetis_direction3[i] = 150 + thetis_direction3[i]
			if thetis_direction3[i] < 140 :
				thetis_direction3[i] = 360 - thetis_direction3[i]

		'''
		for i in range(len(thetis_direction3)):
			if i == 251 :
				thetis_direction3[i] = 150 + thetis_direction3[i]
			if i == 441 :
				thetis_direction3[i] = 150 + thetis_direction3[i]
			if thetis_direction3[i] < 140 :
				thetis_direction3[i] = 360 - thetis_direction3[i]
		b_max = 0
		index_b=0
		for i in range(3):
			for j in range(151*i,151*(i+1)):
				a = abs(velocity_magnitude[j]-thetis_velocity3[j])
				
				b = a/velocity_magnitude[j]
				print(a,b,j)
				if b > b_max:
					b_max = b
					print('b_max',b_max)
					index_b = j
print(b_max,index_b,velocity_magnitude[index_b],thetis_velocity3[index_b])



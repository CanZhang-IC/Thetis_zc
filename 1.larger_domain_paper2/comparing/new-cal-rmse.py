
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
	'5min-16cores-220512-497min'
      ]

for filename in names:
	stationnames=['B5']
	for stationname in stationnames:
		print(filename)
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

		det_file = "../../../outputs/"+filename+"/diagnostic_detectors.hdf5"
		df = h5py.File(det_file, 'r+')
		xvelocity=[]
		yvelocity=[]
		for name, data in df.items():
			if name == stationname:
				xvelocity.append(data[:,1])
				yvelocity.append(data[:,2])
		
		thetis_xvelocity3=[]
		thetis_yvelocity3=[]
		if filename == 'vsotf':
			for i in range(147,198):
				thetis_xvelocity3.append(xvelocity[0][i])
				thetis_yvelocity3.append(yvelocity[0][i])
			for i in range(298,349):
				thetis_xvelocity3.append(xvelocity[0][i])
				thetis_yvelocity3.append(yvelocity[0][i])
			for i in range(482,533):
				thetis_xvelocity3.append(xvelocity[0][i])
				thetis_yvelocity3.append(yvelocity[0][i])
		else:
			# for i in range(2040,2341,2):
			# 	thetis_xvelocity3.append(np.mean(xvelocity[0][i:i+2]))
			# 	thetis_yvelocity3.append(np.mean(yvelocity[0][i:i+2]))
			# for i in range(2952,3253,2):
			# 	thetis_xvelocity3.append(np.mean(xvelocity[0][i:i+2]))
			# 	thetis_yvelocity3.append(np.mean(yvelocity[0][i:i+2]))
			# for i in range(4056,4357,2):
			# 	thetis_xvelocity3.append(np.mean(xvelocity[0][i:i+2]))
			# 	thetis_yvelocity3.append(np.mean(yvelocity[0][i:i+2]))
			for i in range(2040,2341,6):
				thetis_xvelocity3.append(np.mean(xvelocity[0][i-3:i+3]))
				thetis_yvelocity3.append(np.mean(yvelocity[0][i-3:i+3]))
			for i in range(2952,3253,6):
				thetis_xvelocity3.append(np.mean(xvelocity[0][i-3:i+3]))
				thetis_yvelocity3.append(np.mean(yvelocity[0][i-3:i+3]))
			for i in range(4056,4357,6):
				thetis_xvelocity3.append(np.mean(xvelocity[0][i-3:i+3]))
				thetis_yvelocity3.append(np.mean(yvelocity[0][i-3:i+3]))

		thetis_velocity=[]
		thetis_direction=[]
		for i in range(len(thetis_xvelocity3)):
			thetis_velocity.append(sqrt(thetis_xvelocity3[i]**2+thetis_yvelocity3[i]**2))
			one_direction = np.arctan2(thetis_xvelocity3[i],thetis_yvelocity3[i])* 180 / np.pi
			if one_direction < 0:
				thetis_direction.append(one_direction + 360)
			else:
				thetis_direction.append(one_direction)

		if filename == 'vsotf':
			###RMSE###
			thetis_rmse_velocity=[]
			thetis_rmse_direction=[]
			#neap,meso,spring
			for i in range(3):
				thetis_rmse_velocity.append(mean_squared_error(np.array(velocity_magnitude[151*i:151*(i+1):3]),np.array(thetis_velocity[51*i:51*(i+1)]), squared=False))
				#angle_rmse: calculate difference first, then RMSE
				direction_difference = (np.array(velocity_direction[151*i:151*(i+1):3]) - np.array(thetis_direction[51*i:51*(i+1)])+180) % 360 -180
				thetis_rmse_direction.append(np.sqrt(np.sum(direction_difference**2)/direction_difference.shape[0]))
			#total
			thetis_rmse_velocity.append(0)#mean_squared_error(np.array(velocity_magnitude),np.array(thetis_velocity), squared=False))
			#angle_rmse: calculate difference first, then RMSE
			#direction_difference = (np.array(velocity_direction) - np.array(thetis_direction)+180) % 360 -180
			thetis_rmse_direction.append(0)#np.sqrt(np.sum(direction_difference**2)/direction_difference.shape[0]))

			print('The rmse of Thetis-{4} simulated velocity  during neap tide is {0:.3f}, during meso tide is {1:.3f}, during spring tide is {2:.3f} and totally is {3:.3f}.'.
				format(thetis_rmse_velocity[0],thetis_rmse_velocity[1],thetis_rmse_velocity[2],thetis_rmse_velocity[3],filename))
			print('The rmse of Thetis-{4} simulated direction during neap tide is {0:.3f}, during meso tide is {1:.3f}, during spring tide is {2:.3f} and totally is {3:.3f}.'.
				format(thetis_rmse_direction[0],thetis_rmse_direction[1],thetis_rmse_direction[2],thetis_rmse_direction[3],filename))
				
			###R^2###
			thetis_r2_velocity=[]
			thetis_r2_direction=[]
			for i in range(3):
				thetis_r2_velocity.append(r2_score(np.array(velocity_magnitude[151*i:151*(i+1):3]),np.array(thetis_velocity[51*i:51*(i+1)])))
			#all
			thetis_r2_velocity.append(0)#r2_score(np.array(velocity_magnitude),np.array(thetis_velocity)))

			print('The r2 of Thetis-{4} simulated velocity  during neap tide is {0:.3f}, during meso tide is {1:.3f}, during spring tide is {2:.3f} and totally is {3:.3f}.'.
				format(thetis_r2_velocity[0],thetis_r2_velocity[1],thetis_r2_velocity[2],thetis_r2_velocity[3],filename))
		else:
			###RMSE###
			thetis_rmse_velocity=[]
			thetis_rmse_direction=[]
			#neap,meso,spring
			for i in range(3):
				thetis_rmse_velocity.append(mean_squared_error(np.array(velocity_magnitude[151*i:151*(i+1):3]),np.array(thetis_velocity[51*i:51*(i+1)]), squared=False))
				#angle_rmse: calculate difference first, then RMSE
				direction_difference = (np.array(velocity_direction[151*i:151*(i+1):3]) - np.array(thetis_direction[51*i:51*(i+1)])+180) % 360 -180
				thetis_rmse_direction.append(np.sqrt(np.sum(direction_difference**2)/direction_difference.shape[0]))
			#total
			thetis_rmse_velocity.append(0)#mean_squared_error(np.array(velocity_magnitude),np.array(thetis_velocity), squared=False))
			#angle_rmse: calculate difference first, then RMSE
			#direction_difference = (np.array(velocity_direction) - np.array(thetis_direction)+180) % 360 -180
			thetis_rmse_direction.append(0)#np.sqrt(np.sum(direction_difference**2)/direction_difference.shape[0]))

			print('The rmse of Thetis-{4} simulated velocity  during neap tide is {0:.3f}, during meso tide is {1:.3f}, during spring tide is {2:.3f} and totally is {3:.3f}.'.
				format(thetis_rmse_velocity[0],thetis_rmse_velocity[1],thetis_rmse_velocity[2],thetis_rmse_velocity[3],filename))
			print('The rmse of Thetis-{4} simulated direction during neap tide is {0:.3f}, during meso tide is {1:.3f}, during spring tide is {2:.3f} and totally is {3:.3f}.'.
				format(thetis_rmse_direction[0],thetis_rmse_direction[1],thetis_rmse_direction[2],thetis_rmse_direction[3],filename))
				
			###R^2###
			thetis_r2_velocity=[]
			thetis_r2_direction=[]
			for i in range(3):
				thetis_r2_velocity.append(r2_score(np.array(velocity_magnitude[151*i:151*(i+1):3]),np.array(thetis_velocity[51*i:51*(i+1)])))
				#thetis_r2_direction.append(r2_score(np.array(velocity_direction[151*i:151*(i+1):3])/360,np.array(thetis_direction[151*i:151*(i+1):3])/360))
			#all
			thetis_r2_velocity.append(0)#r2_score(np.array(velocity_magnitude),np.array(thetis_velocity)))
			#thetis_r2_direction.append(r2_score(np.array(velocity_direction)/360,np.array(thetis_direction)/360))

			print('The r2 of Thetis-{4} simulated velocity  during neap tide is {0:.3f}, during meso tide is {1:.3f}, during spring tide is {2:.3f} and totally is {3:.3f}.'.
				format(thetis_r2_velocity[0],thetis_r2_velocity[1],thetis_r2_velocity[2],thetis_r2_velocity[3],filename))
			#print('The r2 of Thetis-{4} simulated direction during neap tide is {0:.3f}, during meso tide is {1:.3f}, during spring tide is {2:.3f} and totally is {3:.3f}.'.
				#format(thetis_r2_direction[0],thetis_r2_direction[1],thetis_r2_direction[2],thetis_r2_direction[3],filename))

#show()		

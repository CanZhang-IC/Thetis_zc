import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import numpy as np
import scipy.interpolate 

# Load the xlsx file
excel_data = pd.read_excel('result.xlsx',sheet_name=1,usecols=[0,1,4,5,7])
# print(excel_data.keys())
# Read the values of the file in the dataframe
x = np.array(excel_data['BE'][:])
y = np.array(excel_data['BE_e'][:])
z = np.array(excel_data['P'][:])
points = np.array(list(zip(x,y)))

# print(len(points))
X, Y = np.mgrid[5:46:10, 5:46:10]
Z = scipy.interpolate.griddata(points,z,(X, Y),'cubic')
print(Z)


fig = plt.figure()
ax = plt.axes(projection='3d')

# X, Y, Z = axes3d.get_test_data(0.05)
# ax.plot_surface(X, Y, Z, cmap='rainbow')
# cset = ax.contour(X, Y, Z, zdir='z', offset=-100, cmap=cm.coolwarm)
cset = ax.contour(X, Y, Z, 5,zdir='x', offset=-5, cmap=cm.coolwarm)
# cset = ax.contour(X, Y, Z, 5,zdir='y')
 
ax.set_xlabel('X')
ax.set_xlim(5, 45)
ax.set_ylabel('Y')
# ax.set_ylim(-40, 40)
ax.set_zlabel('Z')
ax.set_zlim(15000, 30000)
 
plt.show()
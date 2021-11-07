import h5py
import numpy

def cal_v(a):
    v = a[1]**2 + a[2]**2
    return numpy.sqrt(v)

df = h5py.File('./diagnostic_detectors.hdf5','r+')
for name,data in df.items():
    if name == 'centre points':
        v_max1 = [1,0,0]
        v_max2 = [-1,0,0]
        for velocity in data:
            # print(velocity)
            if cal_v(velocity) > cal_v(v_max1):
                for ii in range(len(velocity)):
                    v_max1[ii] = velocity[ii]
            elif cal_v(velocity) > cal_v(v_max2) and velocity[1]<0:
                for ii in range(len(velocity)):
                    v_max2[ii] = velocity[ii]

yaw1 = numpy.arctan2(v_max1[1],v_max1[2])*180/numpy.pi
yaw2 = numpy.arctan2(v_max2[1],v_max2[2])*180/numpy.pi
print(cal_v(v_max1),cal_v(v_max2))
print(yaw1,yaw2)
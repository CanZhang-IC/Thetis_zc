#!/bin/bash

while true;
do
for i in 0.4
# do
# cp ../../../outputs/6.yaw_environment/Paper3/Zhoushan_mesh/optimisation/backhome-two_effected/intermediate-op-l_y-P_factor_$i-5min_e\&v/diagnostic_controls.hdf5 ./
# cp ly-forward.py $i.py
# mpirun -np 8 python $i.py;
# rm $i.py
# rm diagnostic_controls.hdf5
# done
for i in 1.0
do
cp ../../../outputs/6.yaw_environment/Paper3/Zhoushan_mesh/optimisation/backhome-two_effected/intermediate-op-l_y-P_factor_$i-5min_e\&v3/diagnostic_controls.hdf5 ./
cp ly-forward.py $i.py
mpirun -np 8 python $i.py;
rm $i.py
rm diagnostic_controls.hdf5
done
sleep 8h;
done


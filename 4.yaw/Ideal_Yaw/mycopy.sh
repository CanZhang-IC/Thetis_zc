for i in 40
do
for j in 0 90
do
    cp angle_optimisation.py angle_optimisation$i-yaw$j.py
    mpirun -np 8 python angle_optimisation$i-yaw$j.py
done
done

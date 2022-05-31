for i in 40
do
for j in 30 60
do
    mpirun -np 8 python angle_optimisation$i-yaw$j.py
done
done

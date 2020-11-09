for i in 10 20 40 60 80 100
do
for j in 0 30
do
    mpirun -np 4 python angle_optimisation$i-yaw$j.py
done
done
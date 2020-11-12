for i in 40
do
for j in 0 30
do
    mpirun -np 8 python angle_optimisation$i-yaw$j.py
done
done

for i in 40
do
for j in 0 30
do
    rm angle_optimisation$i-yaw$j.py
done
done
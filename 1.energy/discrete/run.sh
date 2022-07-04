#! /bin/bash
# for i in $( seq 20 10 40)
# do 
#     cp all-op.py $i.py
#     OMP_NUM_THREADS=1 mpirun -np 4 python $i.py &
#     #rm $i.py
# done
# OMP_NUM_THREADS=1 mpirun -np 4 python all-op.py
# OMP_NUM_THREADS=1 mpirun -np 4 python rectangular-op-yaw.py
OMP_NUM_THREADS=1 mpirun -np 4 python staggered-op-yaw.py

# #! /bin/bash

for i in $( seq 10 -1 0 )
do 
    cp cable_op.py $i.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i.py 
    rm $i.py
    cp cable_forward.py $i.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i.py 
    rm $i.py
done

# 12848
#! /bin/bash
for i in $( seq 10 40 100)
do 
    cp cable_op.py $i.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i.py &
    #rm $i.py
done
# OMP_NUM_THREADS=1 mpirun -np 4 python cable_op.py

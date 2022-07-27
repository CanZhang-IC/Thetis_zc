#! /bin/bash
for i in $( seq 5 5 45)
do 
    cp intermediate.py $i.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i.py 
    rm $i.py
done

# OMP_NUM_THREADS=1 mpirun -np 4 python intermediate.py
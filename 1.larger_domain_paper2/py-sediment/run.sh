#! /bin/bash
for i in $( seq 40 20 40)
do 
    cp run-thetis-sediment.py $i.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i.py 
    rm $i.py
done
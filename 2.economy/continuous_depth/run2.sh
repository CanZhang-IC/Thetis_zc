#! /bin/bash
for i in $( seq 5 5 50)
do 
    cp intermediate_depth.py $i.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i.py
    rm $i.py
done


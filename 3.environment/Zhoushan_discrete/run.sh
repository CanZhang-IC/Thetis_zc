# #! /bin/bash

for i in $( seq 2 2 8 )
do
    cp op-y.py $i.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i.py 
    rm $i.py

    cp forward-y.py $i.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i.py 
    rm $i.py
done

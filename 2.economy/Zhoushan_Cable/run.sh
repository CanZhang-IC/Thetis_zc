# #! /bin/bash

for i in $( seq 0 2 0 )
do
    cp f-e_op.py $i.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i.py 
    rm $i.py

    cp forward_f-e.py $i.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i.py 
    rm $i.py
done

# #! /bin/bash

for i in $( seq 2 2 2 )
do
    cp y_op.py $i.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i.py 
    rm $i.py

    cp forward_y.py $i.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i.py 
    rm $i.py
done

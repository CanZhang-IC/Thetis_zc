# #! /bin/bash

for i in $( seq 2 2 2 )
do
    # cp f-e_op.py $i.py
    # OMP_NUM_THREADS=1 mpirun -np 4 python $i.py 
    # rm $i.py

    cp forward_f-e.py $i.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i.py 
    rm $i.py
done

# #! /bin/bash

for i in $( seq 2 1 2 )
do
    # cp op-1-l.py $i.py
    # OMP_NUM_THREADS=1 mpirun -np 4 python $i.py 
    # rm $i.py

    cp op-2-l_y.py $i.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i.py 
    rm $i.py

    cp op-3-forward.py $i.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i.py 
    rm $i.py
done
# python Hybrid_Code.py
#25294 10129
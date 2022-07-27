#! /bin/bash
# for i in $( seq 7 2 15)
# do 
#     cp neap.py $i.py
#     OMP_NUM_THREADS=1 mpirun -np 4 python $i.py 
#     rm $i.py
# done

OMP_NUM_THREADS=1 mpirun -np 4 python neap.py
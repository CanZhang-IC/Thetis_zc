#! /bin/bash
# for j in $( seq 5 10 50)
# do
# for i in $( seq 5 10 50)
# do 
#     cp forward-behind.py $i-$j.py
#     OMP_NUM_THREADS=1 mpirun -np 4 python $i-$j.py 
#     rm $i-$j.py
# done
# done
OMP_NUM_THREADS=1 mpirun -np 4 python forward-behind.py
#! /bin/bash
for i in $( seq 0 2 8 )
do
for j in $( seq 0 2 8 )
do 
    # cp f-e-op.py $i-$j.py
    # OMP_NUM_THREADS=1 mpirun -np 4 python $i-$j.py 
    # rm $i-$j.py

    cp forward.py $i-$j.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i-$j.py 
    rm $i-$j.py
done
done

# OMP_NUM_THREADS=1 mpirun -np 4 python 6-4.py 
# rm 6-4.py

# cp forward.py 6-4.py
# OMP_NUM_THREADS=1 mpirun -np 4 python 6-4.py 
# rm 6-4.py
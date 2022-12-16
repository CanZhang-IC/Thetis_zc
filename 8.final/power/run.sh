#! /bin/bash
for i in $( seq 40 20 40)
do
for j in $( seq 50 10 50)
do 
    cp f-e-op.py $i-$j.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i-$j.py 
    rm $i-$j.py

    # cp forward-nextto.py $i-$j.py
    # OMP_NUM_THREADS=1 mpirun -np 4 python $i-$j.py 
    # rm $i-$j.py
done
done

# OMP_NUM_THREADS=1 mpirun -np 4 python f-e-op.py 
#! /bin/bash
for j in $( seq 45 20 45)
do
for i in $( seq 75 25 75)
do 
    cp intermediate-behind.py $i-$j.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i-$j.py 
    rm $i-$j.py

    cp forward-behind.py $i-$j.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i-$j.py 
    rm $i-$j.py
done
done
# OMP_NUM_THREADS=1 mpirun -np 4 python forward-behind.py

# 
#! /bin/bash
for j in $( seq 5 20 45)
do
for i in $( seq 25 25 75)
do 
    # cp intermediate-behind.py $i-$j.py
    # OMP_NUM_THREADS=1 mpirun -np 4 python $i-$j.py 
    # # rm $i-$j.py

    cp intermediate-forward.py $i-$j.py
    OMP_NUM_THREADS=1 mpirun -np 4 python $i-$j.py &
    # rm $i-$j.py
done
sleep 0.1h
done

for j in $( seq 5 20 45)
do
for i in $( seq 25 25 75)
do 
    rm $i-$j.py
done
done

#! /bin/bash
for i in $( seq 60 20 100)
do 
    cp intermediate.py $i.py
    mpirun -np 4 python $i.py
    rm $i.py
done


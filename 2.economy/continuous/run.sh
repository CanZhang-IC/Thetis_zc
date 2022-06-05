#! /bin/bash
for i in $( seq 10 10 100)
do 
    cp neap.py $i.py
    mpirun -np 8 python $i.py
    rm $i.py
done
for i in $( seq 10 10 100)
do 
    cp intermediate.py $i.py
    mpirun -np 8 python $i.py
    rm $i.py
done
for i in $( seq 10 10 100)
do 
    cp spring.py $i.py
    mpirun -np 8 python $i.py
    rm $i.py
done
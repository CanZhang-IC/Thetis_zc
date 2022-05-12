#! /bin/bash
for i in $( seq 10 10 100)
do 
    cp intermediate.py intermediate$i.py
    mpirun -np 8 python intermediate$i.py
    rm intermediate$i.py
done
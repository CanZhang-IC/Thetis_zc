#! /bin/bash

for i in $( seq 0 10 360)
do
    cp test7-21.py test7-21$i.py
    mpirun -np 8 python test7-21$i.py
    rm test7-21$i.py
done


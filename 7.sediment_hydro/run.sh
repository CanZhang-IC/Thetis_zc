#! /bin/bash
#$( seq 15 1 17)
for i in  $( seq 5 2 15)
do
    cp ideal_case_run.py $i.py
    mpirun -np 4 python $i.py
    rm $i.py
done

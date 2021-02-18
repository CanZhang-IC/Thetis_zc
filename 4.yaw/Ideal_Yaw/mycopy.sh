#! /bin/bash
for i in $( seq -90 10 90)
do
    cp constant_forward_run.py constant_forward_run$i.py
    mpirun -np 8 python constant_forward_run$i.py
    rm constant_forward_run$i.py
done

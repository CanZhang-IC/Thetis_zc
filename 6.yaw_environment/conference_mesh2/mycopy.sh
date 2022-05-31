#! /bin/bash
for i in $( seq 0.0 0.1 1)
do
    cp forward.py forward$i.py
    mpirun -np 8 python forward$i.py
    rm forward$i.py
done


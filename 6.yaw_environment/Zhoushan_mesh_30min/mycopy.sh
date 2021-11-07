#! /bin/bash
for i in $( seq 0.0 0.1 0.9)
do
    cp intermediate-forward-l_y.py intermediate-forward-l_y$i.py
    mpirun -np 8 python intermediate-forward-l_y$i.py
    rm intermediate-forward-l_y$i.py
done


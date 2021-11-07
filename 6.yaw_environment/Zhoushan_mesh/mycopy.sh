#! /bin/bash
for i in $( seq 0.5 0.1 0.8)
# for i in 0.5 0.1 1
do
    cp intermediate-forward-l_y.py intermediate-forward-l_y$i.py
    mpirun -np 8 python intermediate-forward-l_y$i.py
    rm intermediate-forward-l_y$i.py
done


#! /bin/bash

for i in $( seq 0 10 40)
do
    cp closed_boundary_upwind.py closed_boundary_upwind$i.py
    mpirun -np 8 python closed_boundary_upwind$i.py
    rm closed_boundary_upwind$i.py
done

# shutdown -h now


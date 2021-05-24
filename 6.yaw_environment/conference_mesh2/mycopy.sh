#! /bin/bash
for i in $( seq 0 0.1 1)
do
    cp closed_boundary_upwind.py closed_boundary_upwind$i.py
    mpirun -np 8 python closed_boundary_upwind$i.py
    rm closed_boundary_upwind$i.py
done
# for i in $( seq 0 10 80)
# do
# for j in $( seq 20 10 100)
# do
#     cp closed_boundary_upwind.py closed_boundary_upwind$i$j.py
#     mpirun -np 8 python closed_boundary_upwind$i$j.py
#     rm closed_boundary_upwind$i$j.py
# done
# done

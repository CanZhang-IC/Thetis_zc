#! /bin/bash
for i in $( seq -90 5 90)
# for i in  -45 -30 -25 25 30 45
# do
#     cp closed_boundary.py closed_boundary$i.py
#     mpirun -np 8 python closed_boundary$i.py
#     rm closed_boundary$i.py
# done
# for i in  -45 -35 -30 30 35 45
do
    cp closed_boundary_upwind.py closed_boundary_upwind$i.py
    mpirun -np 8 python closed_boundary_upwind$i.py
    rm closed_boundary_upwind$i.py
done

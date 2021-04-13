#! /bin/bash

mpirun -np 8 python op-closed_boundary_upwind.py

shutdown -h +10 


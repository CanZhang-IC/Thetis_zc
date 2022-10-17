#! /bin/bash

for i in $( seq 10 10 30)
do
    cp one_turbine.py one_turbine$i.py
    OMP_NUM_THREADS=1 mpirun -np 8 python one_turbine$i.py
    rm one_turbine$i.py
    
    # cp two_turbines.py two_turbines$i.py
    # mpirun -np 8 python two_turbines$i.py
    # rm two_turbines$i.py
done


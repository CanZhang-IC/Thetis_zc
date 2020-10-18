for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
do
    mpirun -np 8 python eachoutput_run$i.py
done

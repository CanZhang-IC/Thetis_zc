for i in $( seq 0 10 180)
do
for j in $( seq 20 10 20)
do
    cp forward_test.py forward_test$i$j.py
    mpirun -np 8 python forward_test$i$j.py
    rm forward_test$i$j.py
done
done
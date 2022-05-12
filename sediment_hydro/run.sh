#! /bin/bash
#$( seq 15 1 17)
for i in  $( seq 4 4 40)
do
    cp ideal_case_run.py $i.py
    python $i.py
    rm $i.py
done

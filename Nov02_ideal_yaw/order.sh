for i in 0 30 45 60 90 120 135 150
do
    nohup python angle_optimisation$i.py > logs/$i.out 2>&1 &
done

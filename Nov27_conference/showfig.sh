while (true)
do

cp ../../outputs/conference/onepoint_optimisation/diagnostic_controls.hdf5 ./
cp ../../outputs/conference/onepoint_optimisation/diagnostic_functional.hdf5 ./
python ../read.py
sleep 1800
done
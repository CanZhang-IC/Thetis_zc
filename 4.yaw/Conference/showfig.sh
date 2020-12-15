while (true)
do

cp ../../outputs/conference/test1-optimisation/diagnostic_controls.hdf5 ./
cp ../../outputs/conference/test1-optimisation/diagnostic_functional.hdf5 ./
python ../read.py
sleep 1800
done
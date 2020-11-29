while (true)
do

cp ../../outputs/conference/test2-optimisation-all/diagnostic_controls.hdf5 ./
cp ../../outputs/conference/test2-optimisation-all/diagnostic_functional.hdf5 ./
python ../read.py
sleep 10
done
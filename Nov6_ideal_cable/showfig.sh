while (true)
do
cp ../../outputs/ideal_cable/only_cable/diagnostic_controls.hdf5 ./
cp ../../outputs/ideal_cable/only_cable/diagnostic_functional.hdf5 ./
python ../prepare_cable/Hybrid_Code.py
python ../read.py
sleep 30
done
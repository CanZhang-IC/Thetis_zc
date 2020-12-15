while (true)
do

cp ../../../outputs/3.environment/op_effectonly/diagnostic_controls.hdf5 ./
cp ../../../outputs/3.environment/op_effectonly/diagnostic_functional.hdf5 ./
python ../../read.py
sleep 30
done
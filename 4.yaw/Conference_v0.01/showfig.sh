while (true)
do

cp ../../../outputs/4.yaw/Conference2/oneturbineoptimisation/diagnostic_controls.hdf5 ./
cp ../../../outputs/4.yaw/Conference2/oneturbineoptimisation/diagnostic_functional.hdf5 ./
python ../../read.py
sleep 1800
done
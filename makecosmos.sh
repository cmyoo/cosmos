#!/bin/sh

if [ $# -gt 0 ]; then
DIR=$1
fi

if [ ${DIR} = "adiabatic_spherical" ]; then

if [ -d ${DIR} ]; then
echo "Directory '${DIR}' already exists.  Please rename it."
else
mkdir -p ${DIR}
cp source/makefile ${DIR}/
cp source/*.h ${DIR}/
cp source/*.cpp ${DIR}/
cp sample/${DIR}/*.cpp ${DIR}/
cp sample/${DIR}/*.d ${DIR}/
cp sample/${DIR}/*.dat ${DIR}/
fi

elif [ ${DIR} = "sample_pert" ]; then

if [ -d ${DIR} ]; then
echo "Directory '${DIR}' already exists.  Please rename it."
else
mkdir -p ${DIR}
cp source/makefile ${DIR}/
cp source/*.h ${DIR}/
cp source/*.cpp ${DIR}/
cp sample/${DIR}/*.cpp ${DIR}/
cp sample/${DIR}/*.d ${DIR}/
cp sample/${DIR}/*.dat ${DIR}/
fi

elif [ ${DIR} = "scalar_iso" ]; then

if [ -d ${DIR} ]; then
echo "Directory '${DIR}' already exists.  Please rename it."
else
mkdir -p ${DIR}
cp source/makefile ${DIR}/
cp source/*.h ${DIR}/
cp source/*.cpp ${DIR}/
cp sample/${DIR}/*.cpp ${DIR}/
cp sample/${DIR}/*.d ${DIR}/
cp sample/${DIR}/*.dat ${DIR}/
fi

else

echo "Since there is no sample corresponding to '${DIR}', this script tries to gather necessary files."
echo "After this step, appropriate initial data should be provided to start your simulation by COSMOS code."

if [ -d ${DIR} ]; then
echo "Directory '${DIR}' already exists.  Please rename it."
else
mkdir -p ${DIR}
cp source/makefile ${DIR}/
cp source/*.h ${DIR}/
cp source/*.cpp ${DIR}/
cp source/*.d ${DIR}/
fi

fi


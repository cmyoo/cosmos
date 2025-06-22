#!/bin/sh

if [ $# -gt 0 ]; then
DIR=$1
else
echo "Please input the name of your project as the argument."
echo "You can currently choose the following project as a demo."
echo "  flat_simplest"
echo "  adiabatic_spherical"
echo "  sample_pert"
echo "  scalr_iso"
exit
fi

if [ ${DIR} = "adiabatic_spherical" ]; then

if [ -d ${DIR} ]; then
echo "Directory '${DIR}' already exists.  Please rename it."
else
mkdir -p ${DIR}

if [ -d ${DIR} ]; then
echo "Directory '${DIR}' is created successfully."
else
echo "Directory creation error."
exit
fi

cp source/makefile ${DIR}/
cp source/*.h ${DIR}/
cp source/*.cpp ${DIR}/
cp sample/${DIR}/*.cpp ${DIR}/
cp sample/${DIR}/*.d ${DIR}/
cp sample/${DIR}/*.dat ${DIR}/

fi

elif [ ${DIR} = "flat_simplest" ]; then

if [ -d ${DIR} ]; then
echo "Directory '${DIR}' already exists.  Please rename it."
else
mkdir -p ${DIR}

if [ -d ${DIR} ]; then
echo "Directory '${DIR}' is created successfully."
else
echo "Directory creation error."
exit
fi

#cp source/makefile ${DIR}/
cp source/cosmos.h ${DIR}/
cp source/cosmos_ipol.cpp ${DIR}/
cp source/cosmos_fluid.cpp ${DIR}/
cp source/cosmos_output.cpp ${DIR}/
cp source/cosmos_boundary.cpp ${DIR}/
cp sample/${DIR}/cosmos.cpp ${DIR}/
cp sample/${DIR}/cosmos_bssn.cpp ${DIR}/
cp sample/${DIR}/cosmos_initial.cpp ${DIR}/
cp sample/${DIR}/makefile ${DIR}/
cp sample/${DIR}/par_ini.d ${DIR}/

fi

elif [ ${DIR} = "sample_pert" ]; then

if [ -d ${DIR} ]; then
echo "Directory '${DIR}' already exists.  Please rename it."
else
mkdir -p ${DIR}

if [ -d ${DIR} ]; then
echo "Directory '${DIR}' is created successfully."
else
echo "Directory creation error."
exit
fi

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

if [ -d ${DIR} ]; then
echo "Directory '${DIR}' is created successfully."
else
echo "Directory creation error."
exit
fi

cp source/makefile ${DIR}/
cp source/*.h ${DIR}/
cp source/*.cpp ${DIR}/
cp sample/${DIR}/*.cpp ${DIR}/
cp sample/${DIR}/*.d ${DIR}/
cp sample/${DIR}/*.dat ${DIR}/
fi

else

if [ -d ${DIR} ]; then
echo "Directory '${DIR}' already exists.  Please rename it."
else
mkdir -p ${DIR}

if [ -d ${DIR} ]; then
echo "Directory '${DIR}' is created successfully."
else
echo "Directory creation error."
exit
fi

echo " "
echo "Since there is no sample corresponding to '${DIR}', this script tries to gather necessary files."
echo "After this step, appropriate initial data should be provided to start your simulation by COSMOS code."

cp source/makefile ${DIR}/
cp source/*.h ${DIR}/
cp source/*.cpp ${DIR}/
cp source/*.d ${DIR}/
fi

fi


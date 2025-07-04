#!/bin/sh


if [ $# -gt 0 ]; then
PNAME=$1
else

SAMPLES=`ls --color=never -q1 sample | awk '($1!="man"){print $1}'`

echo "Please input the name of your project."
echo "You can currently choose the following project as a demo or create a new own project."
echo "----------------------------"
echo "### current list of demo ###"
for MODEL in ${SAMPLES};
do
echo ${MODEL}
done
echo "############################"
echo "----------------------------"
echo "If you a new user, please try flat_simplest."
echo "project name ? > "
read PNAME

fi


SAMPLES=`ls --color=never -q1 sample | awk '($1!="man"){print $1}'`

for MODEL in ${SAMPLES};
do

    if [ -z "${MODEL}" ]; then
	break
    fi

    if [ ${PNAME} = ${MODEL} ]; then

	echo "Copying project ${PNAME}..."

	if [ ${PNAME} = "flat_simplest" ]; then

	    if [ -d ${PNAME} ]; then
		echo "Directory '${PNAME}' already exists.  Please rename it before starting the project."
	    else
		mkdir -p ${PNAME}

		if [ -d ${PNAME} ]; then
		    echo "Directory '${PNAME}' is created successfully."
		else
		    echo "Directory creation error."
		    break
		fi

		cp source/cosmos.h ${PNAME}/
		cp source/cosmos_ipol.cpp ${PNAME}/
		cp source/cosmos_fluid.cpp ${PNAME}/
		cp source/cosmos_output.cpp ${PNAME}/
		cp source/cosmos_boundary.cpp ${PNAME}/
		cp sample/${PNAME}/cosmos.cpp ${PNAME}/
		cp sample/${PNAME}/cosmos_bssn.cpp ${PNAME}/
		cp sample/${PNAME}/cosmos_initial.cpp ${PNAME}/
		cp sample/${PNAME}/makefile ${PNAME}/
		cp sample/${PNAME}/par_ini.d ${PNAME}/
		cp sample/${PNAME}/show.gpl ${PNAME}/

		MODEL=""
		break
	    fi

	elif [ ${PNAME} = "flat_simplest_scaleup" ]; then

	    if [ -d ${PNAME} ]; then
		echo "Directory '${PNAME}' already exists.  Please rename it before starting the project."
	    else
		mkdir -p ${PNAME}

		if [ -d ${PNAME} ]; then
		    echo "Directory '${PNAME}' is created successfully."
		else
		    echo "Directory creation error."
		    break
		fi

		cp source/cosmos.h ${PNAME}/
		cp source/cosmos_ipol.cpp ${PNAME}/
		cp source/cosmos_fluid.cpp ${PNAME}/
		cp source/cosmos_boundary.cpp ${PNAME}/
		cp sample/${PNAME}/cosmos.cpp ${PNAME}/
		cp sample/${PNAME}/cosmos_bssn.cpp ${PNAME}/
		cp sample/${PNAME}/cosmos_output.cpp ${PNAME}/
		cp sample/${PNAME}/cosmos_initial.cpp ${PNAME}/
		cp sample/${PNAME}/makefile ${PNAME}/
		cp sample/${PNAME}/par_ini.d ${PNAME}/
		cp sample/${PNAME}/show.gpl ${PNAME}/

		MODEL=""
		break
	    fi

	elif [ ${PNAME} = "flat_simplest_fmr" ]; then

	    if [ -d ${PNAME} ]; then
		echo "Directory '${PNAME}' already exists.  Please rename it before starting the project."
	    else
		mkdir -p ${PNAME}

		if [ -d ${PNAME} ]; then
		    echo "Directory '${PNAME}' is created successfully."
		else
		    echo "Directory creation error."
		    break
		fi

		cp source/cosmos.h ${PNAME}/
		cp source/cosmos_ipol.cpp ${PNAME}/
		cp source/cosmos_fluid.cpp ${PNAME}/
		cp source/cosmos_output.cpp ${PNAME}/
		cp source/cosmos_boundary.cpp ${PNAME}/
		cp sample/${PNAME}/cosmos.cpp ${PNAME}/
		cp sample/${PNAME}/cosmos_fmr.cpp ${PNAME}/
		cp sample/${PNAME}/cosmos_bssn.cpp ${PNAME}/
		cp sample/${PNAME}/cosmos_initial.cpp ${PNAME}/
		cp sample/${PNAME}/makefile ${PNAME}/
		cp sample/${PNAME}/par_ini.d ${PNAME}/
		cp sample/${PNAME}/par_fmr.d ${PNAME}/
		cp sample/${PNAME}/show.gpl ${PNAME}/

		MODEL=""
		break
	    fi

	else

	    if [ -d ${PNAME} ]; then
		echo "Directory '${PNAME}' already exists.  Please rename it before starting the project."
	    else
		mkdir -p ${PNAME}
		
		if [ -d ${PNAME} ]; then
		    echo "Directory '${PNAME}' is created successfully."
		else
		    echo "Directory creation error."
		    break
		fi

		cp source/* ${PNAME}/
		cp sample/${PNAME}/* ${PNAME}/
		MODEL=""
		break
	    fi

	fi

    fi

done


if [ -n "${MODEL}" ];then
    
### if there is no corresponding sample, a new project will be created.

    echo "Creating project name ${PNAME}..."


    if [ -d ${PNAME} ]; then
	echo "Directory '${PNAME}' already exists.  Please rename it."
    else
	mkdir -p ${PNAME}

	if [ -d ${PNAME} ]; then
	    echo "Directory '${PNAME}' is created successfully."
	else
	    echo "Directory creation error."
	    break
	fi

	echo " "
	echo "Since there is no sample corresponding to '${PNAME}', this script tries to gather necessary files."
	echo "After this step, appropriate initial data should be provided to start your simulation by COSMOS code."

	cp source/* ${PNAME}/
	break
    fi
    
fi

#!/bin/sh

if [ -e "out_diff0.dat" ];then

CHECK=`cat out_diff0.dat`
if [ $CHECK -eq 1 ];then
echo "Tests passed."
else
echo "Tests failed.  Check the output file : out_diff.dat. "
fi

else
echo "No expected outputs..."
fi

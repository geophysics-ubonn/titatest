#!/bin/bash
# The git branch, commit, and the host machine are used to generate the binary
# names for CRMod, CRTomo, and CutMcK
# Note: This script generates a Makefile.am file needed by automake!

# Additionally, we compile the list of errors only at autogen-time (i.e. in
# this script) and generate a corresponding crerror.h file to be compiled into
# the binaries.

# Run in src/

###########################
## create crerror.h file ##
###########################

n=$(wc -l error.txt|awk '{printf("%d\n",$1/3)}')
echo "counted $n CR exceptions"

echo $(pwd)

outp='crerror.h'
inp='error.txt'

if [ ! -e $inp ];then
	echo "Error text file not existend"
	exit
fi

echo "CHARACTER(256) :: fetxt($n)" > $outp
echo '' >> $outp

for ((i=1;i<n;i++));do

    ci=$(echo $i|awk '{printf("%03d",$1)}')

    grep $ci $inp|awk '{print}' > tmp

# sed is not used to replace the text (-in) because -in is much slower
    sed "s/$ci//g" tmp > tmp2 # removing line numbers 
    sed '$!N;s/\n//g' tmp2 > tmp # removing end lines
    sed '$!N;s/\n//g' tmp > tmp2 # removing end lines (twice, to remove too much blanks here..)
    cdump=$(paste -s tmp2)     # to be sure paste it--

#    echo "error $ci: $cdump"
    echo "   fetxt($i)='$cdump'" >> $outp
done

rm tmp tmp2

#################################################################
## enter git branch and hostname into the Makefile.am template ##
#################################################################

branch=$( git branch | awk '/\*/{print $2}' )

arch=$( uname -n )

if [ ! -z "${CRT_PREFIX}" ]; then
	prefix="${CRT_PREFIX}"
else
	prefix=$branch'_'$arch
	# replace all '-' characters by '_'
	prefix=${prefix//-/_}
fi

cur=$(pwd)
echo "setting PREFIX $prefix in $cur/Makefile.am"

myorig=./Makefile.orig

mybuff=$(grep myBranchName $myorig|wc -l)

if [ -e $myorig ] && [ ! $mybuff -eq 0 ] ;then
	sed "s/myBranchName/$prefix/g" $myorig > Makefile.am
else
    if [ ! -e $myorig ];then
	echo "no $myorig???"
    else
	echo "myBranchName not found $mybuff"
    fi
fi

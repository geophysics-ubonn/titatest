#!/bin/bash


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
    sed '$!N;s/\n//g' tmp > tmp2 # removing end lines (twice, to remove to much blanks here..)
    cdump=$(paste -s tmp2)     # to be sure paste it--

#    echo "error $ci: $cdump"
    echo "   fetxt($i)='$cdump'" >> $outp
done

rm tmp tmp2

branch=$( git branch | awk '/\*/{print $2}' )

arch=$( uname -n )

prefix=$branch'_'$arch
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

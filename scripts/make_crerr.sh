#!/bin/bash


n_lines=`wc -l error.txt|awk '{printf("%d\n",$1/3)}'`
echo "counted $n_lines CR exceptions"
#n_lines=1
let i=0
let n=$n_lines
outp='crerror.h'

echo "CHARACTER(256) :: fetxt($n_lines)" > $outp
echo '' >> $outp
while [ $i -ne $n ];do
    let i=$i+1
    #echo proof $i
    if [ $i -lt 10 ];then
	ub='00'$i
    elif [ $i -lt 100 ];then
	ub='0'$i
    else
	ub=$i
    fi
    awk -v nl=$ub '{if ($1 == nl){print}}' error.txt > tmp

    sed "s/$ub//g" tmp > tmp2 # removing line numbers
    sed '$!N;s/\n//g' tmp2 > tmp # removing end lines
    sed '$!N;s/\n//g' tmp > tmp2 # removing end lines
    bla=`paste -s tmp2`        # to be sure paste it--

    echo "   fetxt($i)='$bla'" >> $outp
done
rm tmp tmp2

#!/bin/bash

ele=`grep -il "'dat.fin'" *.for`

for x in $ele;do
	sed "s/INCLUDE\ 'dat.fin'/USE\ datmod/g" $x > $x'_2'
	mv $x'_2' $x
done

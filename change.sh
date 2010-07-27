#!/bin/bash

ele=`grep -il "'inv.fin'" *.for`

for x in $ele;do
	sed "s/INCLUDE\ 'inv.fin'/USE\ invmod/g" $x > $x'_2'
	vi $x'_2'
	mv $x'_2' $x
done

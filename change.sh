#!/bin/bash

ele=`grep -il "'sigma.fin'" *.for`

for x in $ele;do
	sed "s/INCLUDE\ 'sigma.fin'/USE\ sigmamod/g" $x > $x'_2'
	vi $x'_2'
	mv $x'_2' $x
done

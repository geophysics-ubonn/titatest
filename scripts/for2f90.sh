#!/bin/bash

for x in *.for;do
	y=`basename $x '.for'`.f90
	sed 's/c\ \ /!!!$\ \ /g' $x > $y
	sed -i 's/c\.\./!!!$\.\./g' $y
	sed -i 's/c::/!!!$::/g' $y
done

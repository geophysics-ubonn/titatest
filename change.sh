#!/bin/bash

#if [ -n "$1" ];then
#	fnm=$1;
#else
#	echo "usage $0 <finname> for substitution";
#	exit;
#fi
ele=`grep -il "'model.fin'" *.for`;

for x in $ele;do
	sed "s/INCLUDE\ 'model.fin'/USE\ modelmod/g" $x > $x'_2'
	vi $x'_2'
	mv $x'_2' $x
done

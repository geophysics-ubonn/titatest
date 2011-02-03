#!/bin/bash

fin="$1"
fout="$2"
liste=`grep -l $fin *.f90`
for x in $liste;do
	sed -i "s/$fin/$fout/g" $x
done

#!/bin/bash

liste=`git branch|sed 's/\*//g' `
echo $liste
for x in $liste;do
	git checkout $x
	git pull
done

#!/bin/bash

modes="Exp Gau Sph"

cur=`pwd`

sgsim="../sgsim"
skel="skel"

for mode in $modes;do
	if [ -d $mode ];then
		echo $mode - directory exists
	else
		mkdir $mode
	fi
	cd $cur/$mode
	sim=`find ./$sgsim -name '*'$mode'.model' -print|sort`
	echo $mode $sim
	for x in $sim;do
		y=`basename "$x" '.model'`
		echo working on $y
		if [ -d "$y" ];then
			echo "rebuild working dir"
			rm -fR $y
		fi
		mkdir $y
		cp -R $cur/$skel/* $y
		cp $x $y/rho/rho.dat
		cd $y/exe
		CRMod_`uname -n` >& crmod_$y.out
		cd $cur
	done
done

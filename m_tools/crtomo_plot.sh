#!/bin/bash

tools="plotCRTmod_batch.m CRTomo epstopdf matlab mogrify";

echo "checking for tools"
forcecheck(){
	Ret_val=`which $1`
	if [ -z "$Ret_val" ];then
		echo "$1 not in search path.. searching $HOME..."
		Ret_val=`find $HOME -name "$1" -print`
		if [ -z "$arg" ];then
			exit
		fi
	fi
	return
}

i=1;
for x in $tools;do
	forcecheck $x
	case "$i" in
		1)
			plotmod=$Ret_val
		;;
		2)
			crt=$Ret_val
		;;
		3)
			epsp=$Ret_val
		;;
		4)
			mtlb=$Ret_val
		;;
		5)
			mogri=$Ret_val
		;;
		*)
			echo "not possible"
			exit
		;;
	esac
	i=$[$i+1];
done

wdir=`pwd|xargs basename`

if [ -n $1 ];then
    wdir=$1
fi

cp $crt crtomo.$wdir

./crtomo.$wdir >& $wdir'.crtrun'

plot_cur_crtomo
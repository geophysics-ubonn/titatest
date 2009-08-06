#!/bin/bash

tools="plotCRTmod_batch.m CRTomo epstopdf matlab mogrify";

forcecheck(){
	echo "checking for $1"
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

rm -f *.eps *.pdf

awk '{if (NR==2){print}}' crtomo.cfg > tmp.meshname

$crt >& $1'.crtrun'
$mtlb < $plotmod >& /dev/null

tmp=`ls *.eps`
myeps=$1'_'$tmp
mv $tmp $myeps
mypdf=`basename $myeps '.eps'`'.pdf'
echo "making image conversion..."
$epsp $myeps -o=$mypdf >& /dev/null
$mogri -format jpg $myeps
echo "done $1"

#pdfcr=`which pdfcrop`
#forcecheck $pdfcr
#$pdfcr --margins 1 --hires $mypdf cropped.pdf >& /dev/null

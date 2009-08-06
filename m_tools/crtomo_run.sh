#!/bin/bash
# 
# Beispiel BASH-script fuer CRTomo lauefe mit unterschiedlichen Skalenlaengen fuer die stochastische regularisierung
# mit auto plot funktion etc.
# Roland Martin 2009
#

crt=`which crtomo_plot.sh`

if [ -z "$crt" ];then
	crt=`which CRTomo`;
	if [ -z "$crt" ];then
		echo "please include either crtomo_plot.sh or CRTomo into your search path.."
		exit
	fi
fi

cur=`pwd`

corrl='1 0.8 0.6 0.5 0.4 0.3 0.2 0.1 0.09 0.08 0.07 0.06 0.05 0.04 0.03 0.02 0.01 0.005 0.001'

for c in $corrl;do
	if [ -d sto_$c ];then
		rm -fR sto_$c;
	fi
	cp -R skel sto_$c
	cd sto_$c;
	awk -v cor=$c '{if (NR==13||NR==14){print cor}else{print}}' crtomo.cfg > tmp.cfg
	mv tmp.cfg crtomo.cfg
	nohup $crt $c &
	cd $cur
done
#!/bin/bash

cur=`pwd`
rho='rho.modl'
ref=$cur/'true'
modes='smo exp gau sph'
fac='100.0'

for x in $modes;do 
	mdir=$cur/$x
	# prepare relative difference plot..
	fnp='pasted_'$x'.modl'
	fnd='diffs_'$x'.modl'
	paste $ref/$rho $mdir/$rho > $fnp 
	awk -v fac=$fac '{if(NR>1){a=(1-($3/$1));printf("%f\t%f\n",log(sqrt(a*a)*fac)/log(10),0.0)}else{print $1}}' $fnp >  $fnd
	l1=`awk -v fac=$fac '{if(NR>1){a=(1-($3/$1));sum=+sqrt(a*a)}} END {printf("%.4f\n",sum/(NR-1)*fac)}' $fnp`
	echo $l1
	clean
	echo "$fnd" > inv.lastmod
	echo "Difference model (1 - true/$x) L_1:$l1" > tmp.fenster
	echo 'log_{10}[%]' > tmp.cbarn
	plot_cur_crtomo
done

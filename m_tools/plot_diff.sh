#!/bin/bash

cur=`pwd`
rho='rho.modl'
vari='inv.variogram'
ref=$cur/'true'
modes='smo exp gau sph'
fac='100.0'



flnd='l1_diff.dat'
flnr='rms.dat'

rm -f $flnd $flnr

for x in $modes;do 
    mdir=$cur/$x
	# prepare relative difference plot..
    fnp='pasted_'$x'.modl'
    fnd='diffs_'$x'.modl'
    paste $ref/$rho $mdir/$rho > $fnp
# get global L1 model difference
    awk -v fac=$fac '{if(NR>1){a=(1-($3/$1));printf("%f\t%f\n",log(sqrt(a*a)*fac)/log(10),0.0)}else{print $1}}' $fnp >  $fnd
    l1=`awk -v fac=$fac '{if(NR>1){a=(1-($3/$1));sum+=sqrt(a*a)}} END {printf("%.2f\n",sum/(NR-1)*fac)}' $fnp`
# paste teh difference to a file
    echo $l1 >> $flnd
# get last RMS and paste to a file
    tail -n 1 $x/inv.stats_it | awk '{print $3}' >> $flnr
    if [ -z $1 ];then
	clean
	echo "$fnd" > inv.lastmod
	echo "Difference model (1 - $x/true) L1=$l1 %" > tmp.fenster
	echo 'log_{10}[%]' > tmp.cbarn
	plot_cur_crtomo >& /dev/null
    fi
done

#!/bin/bash

cur=`pwd`
rho='rho.modl'
vari='inv.variogram'
ref=$cur/'true'
modes='smo exp gau sph'
fac='100.0'



flnd='l1_diff.dat'
flnc='cov.dat'
flnk='korr.dat'
flnr='rms.dat'

rm -f $flnd $flnr $flnc $flnk

for x in $modes;do 
    mdir=$cur/$x
	# prepare relative difference plot..
    fnp='pasted_'$x'.modl'
    fnd='diffs_'$x'.modl'
    paste $ref/$rho $mdir/$rho > $fnp
# get global L1 model difference
    awk -v fac=$fac '{if(NR>1){a=(1-($3/$1));printf("%f\t%f\n",log(sqrt(a*a)*fac)/log(10),0.0)}else{print $1}}' $fnp >  $fnd
    l1=`awk -v fac=$fac '{if(NR>1){a=(1-($3/$1));sum+=sqrt(a*a)}} END {printf("%.1f\n",sum/(NR-1)*fac)}' $fnp`
    means=`awk '{if(NR>1){mx+=$1;my+=$3;n++}}END{print mx/n,my/n}' $fnp`
    mx=`echo $means|awk '{print $1}'`
    my=`echo $means|awk '{print $2}'`
    cov=`awk -v mx=$mx -v my=$my '{if(NR>1){prod+=$1*$3;n++}} END {printf("%.1f\n",(prod-(n*mx*my))/(n-1))}' $fnp`
    korr=`awk -v mx=$mx -v my=$my -v cov=$cov '{if(NR>1){x=$1-mx;y=$3-my;varx+=x*x;vary+=y*y;n++}} END {varx=varx/n;vary=vary/n;printf("%.3f\n",cov/(sqrt(varx)*sqrt(vary)))}' $fnp`
# paste teh difference to a file
    echo $l1 >> $flnd
    echo $cov >> $flnc
    echo $korr >> $flnk
# get last RMS and paste to a file
    tail -n 1 $x/inv.stats_it | awk '{print $3}' >> $flnr
    if [ -z $1 ];then
	clean
	echo "$fnd" > inv.lastmod
	echo "(x=$x) L_1(x,true)=$l1%, Corr(x,true)=$korr" > tmp.fenster
	echo 'log(L_1)[%]' > tmp.cbarn
	echo '-1 3' > tmp.crange
	echo '#' > tmp.noelec
	plot_cur_crtomo >& /dev/null
    fi
done

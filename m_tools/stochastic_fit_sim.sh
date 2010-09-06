#!/bin/bash

echo running $0 with PID $$ at `uname -n`
date

if [ -z "$1" ];then
    var=''
    echo full calculation
else
    var="$1"
    echo only plotting
fi

regud="true smo exp gau sph"
regu=(0 0 15 15 15)

cur=`pwd`

skel=$cur/skel

mode=`basename $cur`

if [ "$mode" == "Exp" ];then
    echo "setting Exponential model"
    rsto=(0 0 0 10 20)
elif [ "$mode" == "Gau" ];then
    echo "setting Gauss model"
    rsto=(1 1 1 11 21)
elif [ "$mode" == "Sph" ];then
    echo "setting Spherical model"
    rsto=(2 2 2 12 22)
else
    exit
fi
for x in *_$mode;do    
    ix=`echo $x|tr '_' ' '|awk '{print $2}'`
    iy=`echo $x|tr '_' ' '|awk '{print $4}'`
# integral scale will be overwritten in this script..
    echo ix:$ix iy:$iy
    cd $x
    crm_skel=crm_skel
    if [ -d "$crm_skel" ];then
	echo $crm_skel exists
    else
	mkdir $crm_skel
	mv config $crm_skel
	mv exe $crm_skel
	mv grid $crm_skel
	mv mod $crm_skel
	mv rho $crm_skel
    fi
    pub=`pwd`/pub
    if [ -d $pub ];then
	echo $pub already there
    else
	mkdir $pub
    fi
    let i=0 # counts the different regus and rstos
    for z in $regud;do
	echo $z $i ${regu[i]} ${rsto[i]}
	if [ -z "$var" ];then
	    if [ -d "$z" ];then
		rm -fR $z;
		rm -fR $pub/$z;
	    fi
	    mkdir $z
	    mkdir $pub/$z
	    cd $z
	    cp -R $cur/$x/$crm_skel/* .
	    cd exe
	    date
	    pwd
	    cp $skel/exe/crtomo.cfg .
	    mybin='CRTomo_'$x'_reg_'$z
	    cp ~/bin/CRTomo_`uname -n` ./$mybin
	    sed -in "s/INVDR/\.\.\/inv/g" crtomo.cfg
	    sed -in "s/RSTO/${rsto[i]}/g" crtomo.cfg
	    sed -in "s/REGU/${regu[i]}/g" crtomo.cfg
	    sed -in "s/IX/$ix/g" crtomo.cfg
	    sed -in "s/IY/$iy/g" crtomo.cfg
	    sed -in "s/ANO/F/g" crtomo.cfg
	    awk '{if(NR==12){print 0}else{print}}' crtomo.cfg > tmp.cfg
	    cp tmp.cfg crtomo.cfg
	    if [ "$z" == "true" ]; then
		awk '{if(NR==8){printf("../rho/rho.dat\n")}else{print}}' crtomo.cfg > tmp.cfg
		awk '{if(NR==15){print 0}else{print}}' tmp.cfg > crtomo.cfg
	    fi
	    if [ ${rsto[i]} -ge 40 ];then
		echo ''|./$mybin >& `echo $mybin|tr 'CRT' 'crt'`'.out'
	    else
		./$mybin >& `echo $mybin|tr 'CRT' 'crt'`'.out'
	    fi
# plot ut the data
	    plt_all_crt $z
	    
	    if [ "$z" == "true" ]; then
# get x-vario fit parameter from true model
		
		cp inv.variogram_x variogram.dat
		vario_fit.sh
		mv variogram_fit.ps variogram_fit_x.ps
		if [ "$mode" == "Exp" ];then
		    ix=`awk '!/#/{print $1}' variofit.dat`
		elif [ $mode == "Gau" ];then 
		    ix=`awk '!/#/{print $2}' variofit.dat`
		elif [ $mode == "Sph" ];then
		    ix=`awk '!/#/{print $3}' variofit.dat`
		fi
# get y-vario fit parameter
		cp inv.variogram_y variogram.dat
		vario_fit.sh
		mv variogram_fit.ps variogram_fit_y.ps
		if [ $mode == "Exp" ];then
		    iy=`awk '!/#/{print $1}' variofit.dat`
		elif [ $mode == "Gau" ];then 
		    iy=`awk '!/#/{print $2}' variofit.dat`
		elif [ $mode == "Sph" ];then
		    iy=`awk '!/#/{print $3}' variofit.dat`
		fi
		echo NEW ix:$ix iy:$iy

		cp inv.variogram variogram.dat
		vario_fit.sh
		mv variogram_fit.ps variogram_fit_r.ps
	    fi
	else
	    cd $z/exe
	fi
	cp inv.variogram* $pub/$z
	cp rho??_mag.pdf $pub/$z/model.pdf
	cp inv.stats_it $pub/$z
	cp inv_stat*.pdf $pub/$z
	cp `cat inv.lastmod|sed 's/mag/modl/g'` $pub/$z/rho.modl
        if [ "$z" == "true" ];then
            cp -R ../grid/ $pub
            sed 's/\.\.\/grid/grid/g' crtomo.cfg > $pub/crtomo.cfg
        fi
        let i=i+1
        cd $cur/$x
    done
    cd $pub
    plot_diff.sh
    variogram.sh $ix $iy $mode
    cd $cur
done

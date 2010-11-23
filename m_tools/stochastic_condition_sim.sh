#!/bin/bash

echo running $0 with PID $$ at `uname -n`
date

get_grid_info(){
    regug=`awk '/Regular/{print}' $1` # if the grid is a regular one, this is in the file..
    echo "checking for a grid info.." 
    if [ -z "$regug" ];then
	echo "there is no regular grid.."
	exit
    else
	echo "found $regug"
    fi
    NX=`awk '{if($1=="NX"){print $2}}' $1`
    X0=`awk '/X0/{print $2}' $1`
    DX=`awk '/DX/{print $2}' $1`
    NY=`awk '/NY/{print $2}' $1`
    Y0=`awk '/Y0/{print $2}' $1`
    DY=`awk '/DY/{print $2}' $1`
    echo "NX:$NX X0:$X0 DX:$DX" # grid settings, for linear equidistant grids
    echo "NY:$NY Y0:$Y0 DY:$DY"
}

get_condition_data (){
	echo "get condition data.."
    if [ -e "$1" ];then
	echo "$1 exists"
    else
	echo "$1 is no valid file"
    fi
    for x in $xlog;do
	for y in $ylog;do
	    echo "cheking value $x $y"
	    val=`awk -v x=$x -v dx=$dxlog -v y=$y -v dy=$dylog -v m=$mean '{if((sqrt( ($1-x)^2 ) < dx) && (sqrt(($2-y)^2) < dy)) { print $1"  "$2"  "$3"  "$3-m } }' $1`
	    if [ -n "$val" ];then
		echo "found $val"
		echo "$val" >> tmp.dat # remove mean from log values
	    fi
	done
    done
    n=`wc -l tmp.dat|awk '{print $1}'`
    echo "found $n matches"
    echo $n	> tmp.head
    echo '4' >> tmp.head
    echo 'x-coord' >> tmp.head
    echo 'y-coord' >> tmp.head
    echo 'resisitivity (log(Ohm m))' >> tmp.head
    echo 'resisitivity-mean (log(Ohm m))' >> tmp.head
    cat tmp.head tmp.dat > $2   
    rm tmp.head tmp.dat
}

if [ -z "$1" ];then
    var=''
    echo full calculation
else
    var="$1"
    echo only plotting
fi

regud="true smo exp gau sph"
regu=(0 0 15 15 15) # regularization switches within crtomo.cfg
# on place of nx
xlog="2.5 7.25 12" # borehole log x-positions
# may be got from inv.elecpositions..
dxlog=".2" # delta x for the information range
ylog="-6" # borehole y-positions
dylog="6" # dy information range 
mean=2.0  #log mean value for modelling-- calc from true.dat using my_meanvariance (grep MEAN) ?
deca=3  # log decades as (maxi(log)-min(log)) -- from true.dat using my_minmax (grep DLOG) ?

cur=`pwd`

skel=$cur/skel
# prepare sgsim parameters
if [ -e "$skel/grid/ctmck.info" ];then
	echo "found $skel/grid/ctmck.info"
	get_grid_info $skel/grid/ctmck.info
else
	echo "trying to run CutMcK.."
	cd $skel/grid
	CutMcK_`uname -n`
	get_grid_info ctmck.info
	cd $cur
fi

mode=`basename $cur`

if [ "$mode" == "Exp" ];then
    echo "setting Exponential model" 
    rsto=(0 0 0 10 20) # this controls the variogram settings for CRTomo, see crtomo.cfg from git repo
    it=1 # variogram model for SGsim
elif [ "$mode" == "Gau" ];then
    echo "setting Gauss model"
    rsto=(1 1 1 11 21)
    it=3
elif [ "$mode" == "Sph" ];then
    echo "setting Spherical model"
    rsto=(2 2 2 12 22)
    it=2
else
    exit
fi

for x in *_$mode;do    
    ix=`echo $x|tr '_' ' '|awk '{print $2}'`
    iy=`echo $x|tr '_' ' '|awk '{print $4}'`
# correlation length will be overwritten in this script..
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
	    # plt_all_crt $z
	    
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
    #plot_diff.sh
    #variogram.sh $ix $iy $mode
    cd $cur
done

#!/bin/bash

minmax (){
    mmin=`awk '!/#/{
	if (NR == 3){p=$3}
	    for(i = 3; i <= NF; i=i+1) { 
		if ($i <= p) {p = $i} 
		    }} END { 
	print p
	}' $1`
    mmax=`awk '!/#/{
	if (NR == 3){m=$3}
	    for(i = 3; i <= NF; i=i+1) { 
		if ($i >= m) {m = $i} 
		    }} END { 
	print m
	}' $1`
    return
}

plotmod () {

    if [ -z "$1" ];then
	echo no argument
	exit
    fi

    if [ -z "$tg1" ];then
	tg1=tmp1
    fi
#
    echo '#' > $tg1
#    echo 'set size 0.34,0.33' >> $tg1
    if [ $i -gt 1 ];then
	echo 'unset key' >> $tg1
    fi
    echo 'set label 4 "'$i')" at graph .07,.94' >> $tg1
    if [ "$ymax" ];then
	echo "set yrange[$ymin:$ymax]" >> $tg1
    fi
    if [ "$ylab" ];then
	echo 'set ylab offset 2 "'$ylab'"' >> $tg1	
    fi
    if [ "$xlab" ];then
#	echo 'set size .34,.55' >> $tg1
	echo 'set xlab offset 0,.5 "'$xlab'"' >> $tg1	
    fi
    echo 'plot \' >> $tg1
    echo '"'$fln'" u 1:3 '$lw' lc 1 ti "smo",\' >> $tg1
    echo '"'$fln'" u 1:4 '$lw' lc 2 ti "exp",\' >> $tg1
    echo '"'$fln'" u 1:5 '$lw' lc 3 ti "gau",\' >> $tg1
    echo '"'$fln'" u 1:6 '$lw' lc 4 ti "sph"' >> $tg1
    cat $tg1 >> $1
    
}

echo running $0 with PID $$ at `uname -n`
date

cur=`pwd`

modes='Exp Gau Sph'


#
tgr=$cur/'tmp.gnu' # plot file..
# declare font size for the legends
font='"Arial,16"' # general font and its size of the title and axes
# major plots
lw='w lp lt 1 lw 3'   # linetype and width of the major plots
lfw='"Times,12"' # legend font
lfw=$font
lfws='spacing 0.8 samplen 0.4'
pw='w p pt 1'
lpw='w l lt 6 lw 2.5'
psw='pointsize 0.5' # pointsize



# MULTIPLOTS!!
#
# set derived gnuplot variables
term="postscript enhanced color font $font"
#term="pngcairo enhanced color crop font $font"
#key="inside bot horizontal font $lfw $lfws"
key="at graph .4,.9 Right vert font $font $lfws width .3 height .2"

# general settings
echo 'set st da l' > $tgr
#echo 'set grid' >> $tgr
echo "set term $term" >> $tgr
#
echo 'set xtics nomirror' >> $tgr
echo 'set ytics nomirror' >> $tgr
echo "set $psw" >> $tgr
echo "set key $key"  >> $tgr

for mode in $modes;do

    flnplot=tmp.$mode.gnu

    if [ -d $mode ];then
	echo working on $mode directory
    else
	echo $mode does not exists.. exiting $0
	exit
    fi

    cd $cur/$mode
    cur2=`pwd`

# set file names containing all the data
    fln_glob=$cur/$mode'_global_l1.dat'
    fln_cov=$cur/$mode'_cov.dat'
    fln_korr=$cur/$mode'_korr.dat'
    fln_rms=$cur/$mode'_rms.dat'
    fln_vr=$cur/$mode'_variogram_l1_r.dat'
    fln_vh=$cur/$mode'_variogram_l1_h.dat'
    fln_vv=$cur/$mode'_variogram_l1_v.dat'
    
    if [ -z "$1" ];then # collecting data
	echo '# L1 global difference evaluation (1-x/true) %' > $fln_glob
	echo '#   Ix    Iy    smo    exp     gau     sph' >> $fln_glob
	
	echo '# Covariance evaluation' > $fln_cov
	echo '#   Ix    Iy    smo    exp     gau     sph' >> $fln_cov

	echo '# Correlation coefficient' > $fln_korr
	echo '#   Ix    Iy    smo    exp     gau     sph' >> $fln_korr

	echo '# RMS evaluation %' > $fln_rms
	echo '#   Ix    Iy    smo    exp     gau     sph' >> $fln_rms
	
	echo '# L1 variogram difference evaluation (r) (1-x/true) %' > $fln_vr
	echo '#   Ix    Iy    smo    exp     gau     sph' >> $fln_vr
	
	echo '# L1 variogram difference evaluation (hori) (1-x/true) %' > $fln_vh
	echo '#   Ix    Iy    smo    exp     gau     sph' >> $fln_vh
	
	echo '# L1 variogram difference evaluation (vert) (1-x/true) %' > $fln_vv
	echo '#   Ix    Iy    smo    exp     gau     sph' >> $fln_vv
	
	
	for x in Ix*_$mode;do
	    ix=`echo $x|tr '_' ' '|awk '{print $2}'`
	    iy=`echo $x|tr '_' ' '|awk '{print $4}'`
	    
	    cd $x
	    
	    pub=`pwd`/pub
	    
	    if [ -d $pub ];then
		echo $pub already there
	    else
		echo nothing to process..
		    exit
	    fi
	    
	    cd $pub
	    
	    plot_diff.sh 1

	    bla=`cat l1_diff.dat`
	    echo $ix $iy $bla >> $fln_glob
	    bla=`cat rms.dat`
	    echo $ix $iy $bla >> $fln_rms
	    bla=`cat cov.dat`
	    echo $ix $iy $bla >> $fln_cov
	    bla=`cat korr.dat`
	    echo $ix $iy $bla >> $fln_korr
	    echo $ix $iy $bla
	    
	    variogram.sh $ix $iy $mode 1
	    
	    cat variogram_l1_r.dat >> $fln_vr
	    cat variogram_l1_h.dat >> $fln_vh
	    cat variogram_l1_v.dat >> $fln_vv
	    
	    cd $cur2
	done
	cd $cur
    else # check if we already perfomed data collection
	if [ -e $fln_rms ];then
	    echo $fln_rms exists, plotting again
	else
	    echo sorry, but there are no results to plot
	    echo run $0 without argument
	    exit
	fi
    fi
    echo '# gnuplot start' > $flnplot
    let i=1
# get glob maxvals RMS (for yrange):
    fln=$fln_rms
    minmax $fln
    ymax=$mmax
    ymin=$mmin
    ymin=`echo $ymin|awk '{print $1-0.2}'`
# RMS plotting command
    ylab='RMS [%]'
    xlab=''
    plotmod $flnplot 'rms'

#
    let i=$i+1
    ymaxv='5.e2'
    if [ "$1" ];then
	ymaxv=$1
    fi
    echo setting upper maximum for L1 diff plot to $ymaxv %
    fln=$fln_glob
    minmax $fln
    ymax=`echo $mmax|awk -v max=$ymaxv '{if($1>max){print max}else{print $1}}'`
    ymin=$mmin
# global L_1 norm of model
    ylab='L_1(x,true) [%]'
    xlab=''
    plotmod $flnplot 'L_1'

    let i=$i+1
# get glob maxvals Correlation
    fln=$fln_korr
    minmax $fln
    ymax=$mmax
    ymin=$mmin
# Covariance plotting command
    ylab='Corr(x,true)'
    xlab=''
    plotmod $flnplot 'corr'

#
    let i=$i+1
    fln=$fln_vr
    minmax $fln
    ymax=$mmax
    ymin=$mmin
    ylab='Global L_1 [%]'
    xlab='Correlation length /[m]'
    plotmod $flnplot 'glob'

    let i=$i+1
    fln=$fln_vh
    minmax $fln
    ymax=$mmax
    ymin=$mmin
    ylab='Horizontal L_1 [%]'
    xlab='Correlation length /[m]'
    plotmod $flnplot 'hor'

    let i=$i+1
    fln=$fln_vv
    minmax $fln
    ymax=$mmax
    ymin=$mmin
    ylab='Vertical L_1 [%]'
    xlab='Correlation length /[m]'
    plotmod $flnplot 'ver'

# now set up special settings for each mplot
# RMS/L_1/corr
    cat $tgr > tmp1
    echo 'set out "'$mode'_eval.ps"'  >> tmp1
    echo 'set origin 0,0' >> tmp1
    echo 'set multiplot layout 2,3 rowsfirst' >> tmp1
#echo 'set size 0.33,0.3' >> tmp1
    cat tmp1 $flnplot > tmp
    gnuplot < tmp

    cd $cur
done


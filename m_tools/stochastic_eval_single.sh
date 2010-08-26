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
    if [ -z "$tg1" ];then
	tg1=tmp2
    fi
    echo "set out $out" > tmp1
    echo "set label 1 $tit $labscr" > $tg1
    if [ "$ymax" ];then
	echo "set yrange[$ymin:$ymax]" >> $tg1
    fi
    if [ "$ylab" ];then
	echo 'set ylab offset 1.5 "'$ylab'"' >> $tg1	
    fi
    if [ "$logy" ];then
	echo 'set log y' >> $tg1
    fi
    echo 'plot \' >> $tg1
    echo '"'$fln'" u 1:3 '$lw' lc 1 ti "smo",\' >> $tg1
    echo '"'$fln'" u 1:4 '$lw' lc 2 ti "exp",\' >> $tg1
    echo '"'$fln'" u 1:5 '$lw' lc 3 ti "gau",\' >> $tg1
    echo '"'$fln'" u 1:6 '$lw' lc 4 ti "sph"' >> $tg1
    
    cat tmp1 $tgr $tg1 > tmp
    
    gnuplot < tmp
    
}
echo running $0 with PID $$ at `uname -n`
date

cur=`pwd`

modes='Exp Gau Sph'

#
tgr=$cur/'tmp.gnu' # plot file..
# declare font size for the legends
font='"Arial,18"' # general font and its size of the title and axes
# major plots
lw='w lp lt 1 lw 3.5'   # linetype and width of the major plots
lfw='"Times,12"' # legend font
lfw=$font
lfws='spacing 0.8 samplen 0.4'
pw='w p pt 1'
lpw='w l lt 6 lw 2.5'
psw='pointsize 0.5' # pointsize

font='"Arial,18"' # general font and its size of the title and axes
# set derived gnuplot variables
term="postscript enhanced color font $font"
#term="pngcairo enhanced color crop font $font"
#key="inside bot horizontal font $lfw $lfws"
key="at screen 0.98,0.98 Right horizontal font $lfw $lfws width 0.4 height 0.2"

# setting global gnuplot parameters 
echo 'set st da l' > $tgr
echo 'set grid' >> $tgr
echo "set term $term" >> $tgr
echo 'set xtics nomirror' >> $tgr
echo 'set xlab offset 0,0.5 "Integral scale /[m]"' >> $tgr
echo 'set ytics nomirror' >> $tgr
echo 'set ylab offset 1.5 "L_1={/Symbol S}_i|1-x_i/y_i| [%]"' >> $tgr
echo "set $psw" >> $tgr
echo "set key $key"  >> $tgr

echo 'set multiplot' >> $tgr
echo 'set origin 0,0' >> $tgr
echo 'set size 1.0,0.97' >> $tgr
labscr='at screen 0.01,0.97'

for mode in $modes;do

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
# get glob maxvals RMS (for yrange):
    fln=$fln_rms
    minmax $fln
    ymax=$mmax
    ymin=$mmin
    ymin=`echo $ymin|awk '{print $1-0.2}'`
# RMS plotting command
    tg1=$fln'.gnu'
    out='"'$fln'.ps"'
    tit='"RMS ('$mode' sim)"'
    ylab='RMS [%]'
    plotmod

# get glob maxvals Covariance (for yrange):
    fln=$fln_cov
    minmax $fln
    ymax=$mmax
    ymin=$mmin
# Covariance plotting command
    tg1=$fln'.gnu'
    out='"'$fln'.ps"'
    tit='"Cov ('$mode' sim)"'
    ylab='Cov(x,y)=({/Symbol S}_ix_iy_i-n@^{\261}x@^{\261}y)/(n-1)'
    plotmod

# get glob maxvals Correlation
    fln=$fln_korr
    minmax $fln
    ymax=$mmax
    ymin=$mmin
# Covariance plotting command
    tg1=$fln'.gnu'
    out='"'$fln'.ps"'
    tit='"Correlation ('$mode' sim)"'
    ylab='Corr(x,y)=Cov(x,y)/((Var(x))^{1/2}(Var(y))^{1/2})'
    plotmod

#
    ymaxv='5.e2'
    if [ "$1" ];then
	ymaxv=$1
    fi
    echo setting upper maximum for L1 diff plot to $ymaxv %
    minmax $fln_glob
    ymax=`echo $mmax|awk -v max=$ymaxv '{if($1>max){print max}else{print $1}}'`
    ymin=$mmin
#    ymax=''
# global L_1 norm of model
    fln=$fln_glob
    tg1=$fln'.gnu'
    out='"'$fln'.ps"'
    tit='"Global model fit ('$mode' sim)"'
    ylab=''
    plotmod

#

    fln=$fln_vr
    minmax $fln
    ymax=$mmax
    ymin=$mmin
    tg1=$fln'.gnu'
    out='"'$fln'.ps"'
    tit='"Global variogram fit ('$mode' sim)"'
    plotmod


    fln=$fln_vh
    minmax $fln
    ymax=$mmax
    ymin=$mmin
    tg1=$fln'.gnu'
    out='"'$fln'.ps"'
    tit='"Horizontal variogram fit ('$mode' sim)"'

    plotmod


    fln=$fln_vv
    minmax $fln
    ymax=$mmax
    ymin=$mmin
    tg1=$fln'.gnu'
    out='"'$fln'.ps"'
    tit='"Vertical variogram fit ('$mode' sim)"'

    plotmod

    cd $cur
done

echo done and cleaning up
rm -f tmp* *gnu
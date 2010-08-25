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

echo running $0 with PID $$ at `uname -n`
date

cur=`pwd`

mode=`basename $cur`

if [ "$mode" == "Exp" ];then
    echo "setting Exponential model 2"
    mtit='Exponential simulations'
elif [ "$mode" == "Gau" ];then
    echo "setting Gauss model"
    mtit='Gauss simulations'
elif [ "$mode" == "Sph" ];then
    echo "setting Spherical model"
    mtit='Spherical simulations'
else
    exit
fi

# prepare data #

# filenames to store data row wise
fln_glob=$cur/$mode'_global_l1.dat'
fln_rms=$cur/$mode'_rms.dat'
fln_vr=$cur/$mode'_variogram_l1_r.dat'
fln_vh=$cur/$mode'_variogram_l1_h.dat'
fln_vv=$cur/$mode'_variogram_l1_v.dat'
    
if [ -z $1 ];then
    echo '# L1 global difference evaluation (1-x/true) %' > $fln_glob
    echo '#   Ix    Iy    smo    exp     gau     sph' >> $fln_glob
    
    echo '# RMS evaluation %' > $fln_rms
    echo '#   Ix    Iy    smo    exp     gau     sph' >> $fln_rms
    
    echo '# L1 variogram difference evaluation (r) (1-x/true) %' > $fln_vr
    echo '#   Ix    Iy    smo    exp     gau     sph' >> $fln_vr
    
    echo '# L1 variogram difference evaluation (hori) (1-x/true) %' > $fln_vh
    echo '#   Ix    Iy    smo    exp     gau     sph' >> $fln_vh
    
    echo '# L1 variogram difference evaluation (vert) (1-x/true) %' > $fln_vv
    echo '#   Ix    Iy    smo    exp     gau     sph' >> $fln_vv
    
    
    for x in *_$mode;do
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
	echo $ix $iy $bla 
	echo $ix $iy $bla >> $fln_rms
	
	variogram.sh $ix $iy $mode 1
	
	cat variogram_l1_r.dat >> $fln_vr
	cat variogram_l1_h.dat >> $fln_vh
	cat variogram_l1_v.dat >> $fln_vv
	
	cd $cur
    done
else
    if [ -e $fln_rms ];then
	echo $fln_rms exists, plotting again
    else
	echo sorry, but there are no results to plot
	echo run $0 without argument
	exit
    fi
fi
#
# MULTIPLOTS!!
#
#
tgr=$cur/'tmp.gnu' # plot file..
# declare font size for the legends
font='"Arial,16"' # general font and its size of the title and axes
# major plots
lw='w lp lt 1 lw 3.5'   # linetype and width of the major plots
lfw='"Times,12"' # legend font
lfw=$font
lfws='spacing 0.8 samplen 0.4'
pw='w p pt 1'
lpw='w l lt 6 lw 2.5'
psw='pointsize 0.5' # pointsize

# set derived gnuplot variables
term="postscript enhanced color font $font"
#term="pngcairo enhanced color crop font $font"
out=$mode'_sto_eval.ps'
#key="inside bot horizontal font $lfw $lfws"
key="at screen 0.5,0.48 center horizontal font $lfw $lfws width 0.4 height 0.2"

# multi
# setting global gnuplot parameters 
echo 'set st da l' > $tgr
echo 'set grid' >> $tgr
echo "set term $term" >> $tgr
# for global model diff plot plot also the RMS
echo 'set y2tics nomirror' >> $tgr
#echo 'set y2label offset -2 "RMS [%]"' >> $tgr
echo 'set log y2' >> $tgr
#
echo 'set xtics nomirror' >> $tgr
echo 'set ytics nomirror' >> $tgr
#echo 'set log y' >> $tgr
echo 'set ylab offset 1.5 "(1-x/true) [%]"' >> $tgr
echo "set $psw" >> $tgr
echo "set key $key"  >> $tgr

echo 'set out "'$out'"' >> $tgr

echo 'set label 1 "'$mtit'" at screen 0.5,0.98 center' >> $tgr # main plot title
echo 'set label 2 "a) Global model fit" at screen 0.02,0.97' >> $tgr
echo 'set label 3 "b) Global variogram fit" at screen 0.02,0.46' >> $tgr
echo 'set label 4 "c) Horizontal variogram fit" at screen 0.735,0.97' >> $tgr
echo 'set label 5 "d) Vertical variogram fit" at screen 0.76,0.46' >> $tgr
#echo 'set label 6 "RMS [%]" at screen 0.46,0.67' >> $tgr

echo 'set multiplot' >> $tgr
echo 'set size 0.5,0.48' >> $tgr # equal sizes..
#
# setting specific values 
#
# top left, global error

echo 'set origin 0,0.5' >> $tgr
echo 'set size 0.58,0.43' >> $tgr
echo 'plot \' >> $tgr
echo '"'$fln_glob'" u 1:3 '$lw' lc 1 ti "smo",\' >> $tgr
echo '"'$fln_rms'" u 1:3 axes x1y2 '$lpw' lc 1 ti "(rms)",\' >> $tgr
echo '"'$fln_glob'" u 1:4 '$lw' lc 2 ti "exp",\' >> $tgr
echo '"'$fln_rms'" u 1:4 axes x1y2 '$lpw' lc 2 ti "(rms)",\' >> $tgr
echo '"'$fln_glob'" u 1:5 '$lw' lc 3 ti "gau",\' >> $tgr
echo '"'$fln_rms'" u 1:5 axes x1y2 '$lpw' lc 3 ti "(rms)",\' >> $tgr
echo '"'$fln_glob'" u 1:6 '$lw' lc 4 ti "sph",\' >> $tgr
echo '"'$fln_rms'" u 1:6 axes x1y2 '$lpw' lc 4 ti "(rms)"' >> $tgr

    
# variogram plots
echo 'unset y2tics' >> $tgr
echo 'unset y2label' >> $tgr
echo 'unset log y' >> $tgr
echo 'set yrange [*:*]' >> $tgr
echo 'unset ylabel' >> $tgr
echo 'unset key' >> $tgr

# top right horizontal variogram error 
echo 'set size 0.52,0.43' >> $tgr # equal sizes..
echo 'set origin 0.5,0.5' >> $tgr
echo 'plot \' >> $tgr
echo '"'$fln_vh'" u 1:3 '$lw' lc 1 ti "smo",\' >> $tgr
echo '"'$fln_vh'" u 1:4 '$lw' lc 2 ti "exp",\' >> $tgr
echo '"'$fln_vh'" u 1:5 '$lw' lc 3 ti "gau",\' >> $tgr
echo '"'$fln_vh'" u 1:6 '$lw' lc 4 ti "sph"' >> $tgr


# bottom left, global variogram error
echo 'set size 0.52,0.46' >> $tgr # equal sizes..
echo 'set xlab offset 0,0.5 "Integral scale /[m]"' >> $tgr
echo 'set ylab offset 1.5 "(1-x/true) [%]"' >> $tgr
echo 'set origin 0,0' >> $tgr
echo 'plot \' >> $tgr
echo '"'$fln_vr'" u 1:3 '$lw' lc 1 ti "smo",\' >> $tgr
echo '"'$fln_vr'" u 1:4 '$lw' lc 2 ti "exp",\' >> $tgr
echo '"'$fln_vr'" u 1:5 '$lw' lc 3 ti "gau",\' >> $tgr
echo '"'$fln_vr'" u 1:6 '$lw' lc 4 ti "sph"' >> $tgr

# bottom right vertical variogram error
echo 'set origin 0.5,0' >> $tgr
echo 'unset ylabel' >> $tgr
echo 'plot \' >> $tgr
echo '"'$fln_vv'" u 1:3 '$lw' lc 1 ti "smo",\' >> $tgr
echo '"'$fln_vv'" u 1:4 '$lw' lc 2 ti "exp",\' >> $tgr
echo '"'$fln_vv'" u 1:5 '$lw' lc 3 ti "gau",\' >> $tgr
echo '"'$fln_vv'" u 1:6 '$lw' lc 4 ti "sph"' >> $tgr

gnuplot < $tgr
my_pscrop $out
mv `echo $out|sed 's/\"//g'|sed 's/\.ps/\.pdf/g'` ..

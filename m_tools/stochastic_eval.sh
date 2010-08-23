#!/bin/bash

echo running $0 with PID $$ at `uname -n`
date

cur=`pwd`

mode=`basename $cur`

fln_glob=$cur/$mode'_global_l1.dat'

fln_vr=$cur/$mode'_variogram_l1_r.dat'
fln_vh=$cur/$mode'_variogram_l1_h.dat'
fln_vv=$cur/$mode'_variogram_l1_v.dat'

if [ "$mode" == "Exp" ];then
    echo "setting Exponential model 2"
elif [ "$mode" == "Gau" ];then
    echo "setting Gauss model"
elif [ "$mode" == "Sph" ];then
    echo "setting Spherical model"
else
    exit
fi

echo '# L1 global difference evaluation (1-x/true) %' > $fln_glob
echo '#   Ix    Iy    smo    exp     gau     sph' >> $fln_glob

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
    echo $ix $iy $bla 
    echo $ix $iy $bla >> $fln_glob

    variogram.sh $ix $iy $mode 1

    cat variogram_l1_r.dat >> $fln_vr
    cat variogram_l1_h.dat >> $fln_vh
    cat variogram_l1_v.dat >> $fln_vv

    cd $cur
done

# declare font size for the legends
font='14' # general font and its size of the title and axes
# major plots
lw='w l lw 4'   # linetype and width of the major plots
lfw='Times,12' # legend font
lfws='spacing 0.9'
pw='w lp lw 3'
psw='pointsize 0.8' # pointsize

tg1='tmp.gnu'
tg2=$mode'_l1_evaluations'

# setting global gnuplot parameters
echo 'set st da l' > $tg1
echo 'set grid' >> $tg1
echo "set term pos enh col sol $font" >> $tg1
echo 'set xtics nomirror' >> $tg1
echo 'set ytics nomirror' >> $tg1

echo 'set xlab offset 0,0.5 "Integral scale /[m]"' >> $tg1
echo 'set ylab offset 2 "L_1 difference (1-x/true) [%]"' >> $tg1
echo "set key inside bot left Right nobox noreverse $lfws" >> $tg1
echo "set $psw" >> $tg1
echo 'set out "'$tg2'.ps"' >> $tg1
echo 'set multiplot' >> $tg1
echo 'set size 0.52,0.52' >> $tg1 # equal sizes..
#
# setting specific values 
#
# top left, global error
echo 'set origin 0,0.5' >> $tg1
tit="Global model fit ($mode)"
echo 'set tit "'$tit'"' >> $tg1
echo 'plot \' >> $tg1
echo '"'$fln_glob'" u 1:3 '$pw' ti "smo",\' >> $tg1
echo '"'$fln_glob'" u 1:4 '$pw' ti "exp",\' >> $tg1
echo '"'$fln_glob'" u 1:5 '$pw' ti "gau",\' >> $tg1
echo '"'$fln_glob'" u 1:6 '$pw' ti "sph"' >> $tg1

echo 'unset key' >> $tg1
# top right global variogram error 
echo 'set origin 0.5,0.5' >> tmp.gnu
tit="Global variogram fit ($mode)"
echo 'set tit "'$tit'"' >> $tg1
echo 'plot \' >> $tg1
echo '"'$fln_vr'" u 1:3 '$pw' ti "smo",\' >> $tg1
echo '"'$fln_vr'" u 1:4 '$pw' ti "exp",\' >> $tg1
echo '"'$fln_vr'" u 1:5 '$pw' ti "gau",\' >> $tg1
echo '"'$fln_vr'" u 1:6 '$pw' ti "sph"' >> $tg1

# bottom left, horizontal variogram error
echo 'set origin 0,0' >> tmp.gnu
tit="Horizontal variogram fit ($mode)"
echo 'set tit "'$tit'"' >> $tg1
echo 'plot \' >> $tg1
echo '"'$fln_vh'" u 1:3 '$pw' ti "smo",\' >> $tg1
echo '"'$fln_vh'" u 1:4 '$pw' ti "exp",\' >> $tg1
echo '"'$fln_vh'" u 1:5 '$pw' ti "gau",\' >> $tg1
echo '"'$fln_vh'" u 1:6 '$pw' ti "sph"' >> $tg1

# bottom right vertical variogram error
echo 'set origin 0.5,0' >> tmp.gnu
tit="Vertical variogram fit ($mode)"
echo 'set tit "'$tit'"' >> $tg1
echo 'plot \' >> $tg1
echo '"'$fln_vv'" u 1:3 '$pw' ti "smo",\' >> $tg1
echo '"'$fln_vv'" u 1:4 '$pw' ti "exp",\' >> $tg1
echo '"'$fln_vv'" u 1:5 '$pw' ti "gau",\' >> $tg1
echo '"'$fln_vv'" u 1:6 '$pw' ti "sph"' >> $tg1


gnuplot < $tg1
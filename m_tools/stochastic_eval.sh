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
font='"Helvetica" 16' # general font and its size of the title and axes
# major plots
lw='w l lw 4'   # linetype and width of the major plots
lfw='Times,12' # legend font
lfws='spacing 1.1'
psw='pointsize 1.5' # pointsize

tg1='tmp.gnu'

echo 'set st da l' > $tg1
echo 'set grid' >> $tg1
echo "set term pos enh col sol $font" >> $tg1
echo 'set xtics nomirror' >> $tg1
echo 'set ytics nomirror' >> $tg1

echo 'set xlab offset 0,0.5 "Integral scale /[m]"' >> $tg1
echo 'set ylab offset 2 "L1 difference [%]"' >> $tg1
echo 'se key inside right bot nobox noreverse font "'$lfw'"' $lfws >> $tg1
echo "set $psw" >> $tg1

echo 'set out "l1_evaluation.ps"' >> $tg1

echo 'set multiplot' >> $tg1
echo 'set size 0.52,0.5' >> $tg1
# top left global error
echo 'set origin 0,0.5' >> $tg1

tit="$mode L1 global model fit against reference"
echo 'set tit "'$tit'"' >> $tg1

echo 'plot \' >> $tg1
echo '"'$mode'_" u 1:($3/true) '$lw' lc 0 ti "{/Symbol g}(h)",\' >> $tg1
# bottom left horizontal variogram error
echo 'set origin 0,0' >> tmp.gnu
# top right global variogram error 
echo 'set origin 0.5,0.5' >> tmp.gnu
echo 'set size 0.5,0.5' >> tmp.gnu
# bottom right vertical variogram error
echo 'set origin 0.48,0' >> tmp.gnu

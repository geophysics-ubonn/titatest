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

echo '# L1 global difference evaluation' > $fln_glob
echo '#   Ix    Iy    smo    exp     gau     sph' >> $fln_glob

echo '# L1 variogram difference evaluation (r)' > $fln_vr
echo '#   Ix    Iy    true   smo    exp     gau     sph' >> $fln_vr

echo '# L1 variogram difference evaluation (hori)' > $fln_vh
echo '#   Ix    Iy    true   smo    exp     gau     sph' >> $fln_vh

echo '# L1 variogram difference evaluation (vert)' > $fln_vv
echo '#   Ix    Iy    true   smo    exp     gau     sph' >> $fln_vv


for x in *_$mode;do
    ix=`echo $x|tr '_' ' '|awk '{print $2}'`
    iy=`echo $x|tr '_' ' '|awk '{print $4}'`
    echo $ix $iy > tmp.scl
    cd $x
    pub=`pwd`/pub
    if [ -d $pub ];then
	echo $pub already there
    else
	echo nothing to process..
	exit
    fi
    cd $pub
    plot_diff.sh

    paste tmp.scl l1_diff.dat > tmp.dat
    cat tmp.dat >> $fln_glob

    variogram.sh $ix $iy $mode

    cat variogram_l1_r.dat >> $vln_vr
    cat variogram_l1_h.dat >> $vln_vh
    cat variogram_l1_v.dat >> $vln_vv

    cd $cur
done

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
    if [ $mode == "Exp" ];then
	echo 'unset xlab' > $tg1
	echo 'set label 1 "a) Exp sim" at screen 0.02,0.97' >> $tg1
	echo 'set size 1.0,0.33' >> $tg1
	echo 'set origin 0,0.64' >> $tg1
    elif [ $mode == "Gau" ];then
	echo 'unset key' > $tg1
	echo 'set size 1.0,0.31' >> $tg1
	echo 'set origin 0,0.33' >> $tg1
	echo 'set label 1 "b) Gau sim" at screen 0.02,0.64' >> $tg1
    elif [ $mode == "Sph" ];then
	echo 'set label 1 "c) Sph sim" at screen 0.02,0.33' > $tg1
	echo 'set xlab offset 0,0.5 "Integral scale /[m]"' >> $tg1
    fi
    if [ "$ymax" ];then
	echo "set yrange[$ymin:$ymax]" >> $tg1
    fi
    if [ "$ylab" ];then
	echo 'set ylab offset 1.5 "'$ylab'"' >> $tg1	
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
font='"Arial,18"' # general font and its size of the title and axes
# major plots
lw='w lp lt 1 lw 3.5'   # linetype and width of the major plots
lfw='"Times,12"' # legend font
lfw=$font
lfws='spacing 0.8 samplen 0.4'
pw='w p pt 1'
lpw='w l lt 6 lw 2.5'
psw='pointsize 0.5' # pointsize



rms=$cur/tmp.rms.gnu
l1_diff=$cur/tmp.l1_diff.gnu
vario_glo=$cur/tmp.vario_glo.gnu
vario_ver=$cur/tmp.vario_ver.gnu
vario_hor=$cur/tmp.vario_hor.gnu

rm -f $rms $l1_diff $vario_glo $vario_ver $vario_hor

for mode in $modes;do

    if [ -d $mode ];then
	echo working on $mode directory
    else
	echo $mode does not exists.. exiting $0
	exit
    fi

    cd $cur/$mode
    cur2=`pwd`
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
	    
	    cd $cur2
	done
    else
	if [ -e $fln_glob ];then
	    echo $fln_rms exists, plotting again
	else
	    echo "usage: $0 <ymaxval>"
	fi
    fi

    # get maxvals of RMS vals
    fln=$fln_rms
    minmax $fln
    ymax=$mmax
    ymin=$mmin
    ymin=`echo $ymin|awk '{print $1-0.2}'`
#    ymax=''
    echo $fln $ymax $ymin
    ylab='RMS [%]'

    plotmod $rms

# L1 diff
    fln=$fln_glob
    ymaxv='5.e2'
    if [ "$1" ];then
	ymaxv=$1
    fi
    echo setting upper maximum for L1 diff plot to $ymaxv %
    minmax $fln
    ymax=`echo $mmax|awk -v max=$ymaxv '{if($1>max){print max}else{print $1}}'`
    ymin=$mmin
    ylab=''

    plotmod $l1_diff

# global variogram plotting commands
    fln=$fln_vr
    minmax $fln
    ymax=$mmax
    ymin=$mmin
    ylab=''
    plotmod $vario_glo

    fln=$fln_vh
    minmax $fln
    ymax=$mmax
    ymin=$mmin
    ylab=''
    plotmod $vario_hor

    fln=$fln_vv
    minmax $fln
    ymax=$mmax
    ymin=$mmin
    ylab=''
    plotmod $vario_ver


#    rm tmp*

    cd $cur
done
# MULTIPLOTS!!
#
# set derived gnuplot variables
term="postscript enhanced color font $font"
#term="pngcairo enhanced color crop font $font"
#key="inside bot horizontal font $lfw $lfws"
key="at screen 1,.97 Right horizontal font $font $lfws width 0.3 height 0.2"

titpos='at screen 0.51,0.98 center'
# general settings
echo 'set st da l' > $tgr
echo 'set grid' >> $tgr
echo "set term $term" >> $tgr
#
echo 'set xtics nomirror' >> $tgr
echo 'set ytics nomirror' >> $tgr
echo 'set ylab offset 1.5 "(1-x/true) [%]"' >> $tgr
echo "set $psw" >> $tgr
echo "set key $key"  >> $tgr


# now set up special settings for each mplot
# RMS
tit='Inversion RMS'
echo 'set out "rms.ps"'  > tmp1
echo 'set label 2 "'$tit'"' $titpos >> tmp1
echo 'set origin 0,0' >> tmp1
echo 'set multiplot layout 3,1' >> tmp1

cat $tgr tmp1 $rms > tmp
gnuplot < tmp


# L1 diff of global model
tit='Global model L_1-difference'
echo 'set out "l1_diff.ps"'  > tmp1
echo 'set label 2 "'$tit'"' $titpos >> tmp1
echo 'set origin 0,0' >> tmp1
echo 'set multiplot layout 3,1' >> tmp1

cat $tgr tmp1 $l1_diff > tmp
gnuplot < tmp

# Variogram
tit='Total variogram fit (L_1) '
echo 'set out "vario_glo.ps"'  > tmp1
echo 'set origin 0,0' >> tmp1
echo 'set label 2 "'$tit'"' $titpos >> tmp1
echo 'set multiplot layout 3,1' >> tmp1
cat $tgr tmp1 $vario_glo > tmp
gnuplot < tmp


tit='Horizontal variogram fit (L_1) '
echo 'set out "vario_hor.ps"'  > tmp1
echo 'set origin 0,0' >> tmp1
echo 'set label 2 "'$tit'"' $titpos >> tmp1
echo 'set multiplot layout 3,1' >> tmp1
cat $tgr tmp1 $vario_hor > tmp
gnuplot < tmp


tit='Vertical variogram fit (L_1) '
echo 'set out "vario_ver.ps"'  > tmp1
echo 'set origin 0,0' >> tmp1
echo 'set label 2 "'$tit'"' $titpos >> tmp1
echo 'set multiplot layout 3,1' >> tmp1
cat $tgr tmp1 $vario_ver > tmp
gnuplot < tmp

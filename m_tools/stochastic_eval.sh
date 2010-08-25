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

# multi
# setting global gnuplot parameters 
echo 'set st da l' > $tgr
echo 'set grid' >> $tgr
echo "set term $term" >> $tgr
# for global model diff plot plot also the RMS
echo 'set y2tics nomirror' >> $tgr
echo 'set y2label offset -1.5 "RMS [%]"' >> $tgr
echo 'set log y2' >> $tgr
#
echo 'set xtics nomirror' >> $tgr
echo 'set xlab offset 0,0.5 "Integral scale /[m]"' >> $tgr
echo 'set ytics nomirror' >> $tgr
#echo 'set log y' >> $tg1
echo 'set ylab offset 1.5 "L_1 difference (1-x/true) [%]"' >> $tgr
echo "set $psw" >> $tgr
echo "set key $key"  >> $tgr

echo 'set multiplot' >> $tgr
echo 'set origin 0,0' >> $tgr
echo 'set size 1.0,0.97' >> $tgr
labscr='at screen 0.01,0.97'

rms_glo=$cur/tmp.l1_rms_glo
vario_glo=$cur/tmp.vario_glo
vario_ver=$cur/tmp.vario_ver
vario_hor=$cur/tmp.vario_hor
rm -f $rms_glo $vario_glo $vario_ver $vario_hor


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
	if [ -e $fln_rms ];then
	    echo $fln_rms exists, plotting again
	else
	    echo sorry, but there are no results to plot
	    echo run $0 without argument
	    exit
	fi
    fi
    # get glob maxvals l1diff (yrange):
    minmax $fln_rms
    y2max=$mmax
    y2min=$mmin
    y2min=`echo $y2min|awk '{print $1-0.2}'`

    ymaxv=5.e2
    if [ -n $1 ];then
	ymaxv=$1
    fi
    echo setting upper maximum for L1 diff plot to $ymaxv %
    minmax $fln_glob
    ymax=`echo $mmax|awk -v max=$ymaxv '{if($1>max){print max}else{print $1}}'`
    ymin=$mmin
#
    if [ $mode == "Exp" ];then
	echo 'set size 1.0,0.33' > tmp.$mode
	echo 'set label 1 "a) Exp" at screen 0.02,0.97' >> tmp.$mode
	echo 'set origin 0,0.64' >> tmp.$mode
    elif [ $mode == "Gau" ];then
	echo 'unset key' > tmp.$mode
	echo 'set size 1.0,0.31' >> tmp.$mode
	echo 'set origin 0,0.33' >> tmp.$mode
	echo 'set label 1 "b) Gau" at screen 0.02,0.65' >> tmp.$mode
    elif [ $mode == "Sph" ];then
	echo 'set label 1 "c) Sph" at screen 0.02,0.34' > tmp.$mode
	echo 'set xlab offset 0,0.5 "Integral scale /[m]"' >> tmp.$mode
    fi

# global L1 fit and RMS plotting command
    fln=$fln_glob
    tg1=tmp1
    out='"'$fln'.ps"'
    tit='"Global model fit and RMS ('$mode')"'
    
    echo "set out $out" > tmp.gnu
    echo "set label 1 $tit $labscr" >> tmp.gnu
#    echo 'unset key' > $tg1
#    echo 'set key inside left Left font "Arial,14"' $lfws >> $tg1  
    echo "set yrange[$ymin:$ymax]" > $tg1
    echo "set y2range [$y2min:$y2max]" >> $tg1
    echo 'plot \' >> $tg1
    echo '"'$fln_glob'" u 1:3 '$lw' lc 1 ti "smo",\' >> $tg1
    echo '"'$fln_rms'" u 1:3 axes x1y2 '$lpw' lc 1 ti "(rms)",\' >> $tg1
    echo '"'$fln_glob'" u 1:4 '$lw' lc 2 ti "exp",\' >> $tg1
    echo '"'$fln_rms'" u 1:4 axes x1y2 '$lpw' lc 2 ti "(rms)",\' >> $tg1
    echo '"'$fln_glob'" u 1:5 '$lw' lc 3 ti "gau",\' >> $tg1
    echo '"'$fln_rms'" u 1:5 axes x1y2 '$lpw' lc 3 ti "(rms)",\' >> $tg1
    echo '"'$fln_glob'" u 1:6 '$lw' lc 4 ti "sph",\' >> $tg1
    echo '"'$fln_rms'" u 1:6 axes x1y2 '$lpw' lc 4 ti "(rms)"' >> $tg1

    echo 'set size 1.0,0.92' > tmp

    cat $tgr tmp $tg1 >> tmp.gnu
    gnuplot < tmp.gnu
    cat tmp.$mode $tg1 >> $rms_glo

# global variogram plotting commands
    fln=$fln_vr
    tg1=tmp1
    out='"'$fln'.ps"'
    tit='"Global variogram fit ('$mode')"'
    
    echo "set out $out" > tmp.gnu
    echo "set label 1 $tit $labscr" >> tmp.gnu

    echo 'unset y2tics' > $tg1
    echo 'unset y2label' >> $tg1
    echo 'plot \' >> $tg1
    echo '"'$fln'" u 1:3 '$lw' lc 1 ti "smo",\' >> $tg1
    echo '"'$fln'" u 1:4 '$lw' lc 2 ti "exp",\' >> $tg1
    echo '"'$fln'" u 1:5 '$lw' lc 3 ti "gau",\' >> $tg1
    echo '"'$fln'" u 1:6 '$lw' lc 4 ti "sph"' >> $tg1

    cat $tgr $tg1 >> tmp.gnu
    gnuplot < tmp.gnu
    cat tmp.$mode $tg1 >> $vario_glo

    fln=$fln_vh
    tg1=tmp1
    out='"'$fln'.ps"'
    tit='"Horizontal variogram fit ('$mode')"'
    
    echo "set out $out" > tmp.gnu
    echo "set label 1 $tit $labscr" >> tmp.gnu

    echo 'unset y2tics' > $tg1
    echo 'unset y2label' >> $tg1
    echo 'plot \' >> $tg1
    echo '"'$fln'" u 1:3 '$lw' lc 1 ti "smo",\' >> $tg1
    echo '"'$fln'" u 1:4 '$lw' lc 2 ti "exp",\' >> $tg1
    echo '"'$fln'" u 1:5 '$lw' lc 3 ti "gau",\' >> $tg1
    echo '"'$fln'" u 1:6 '$lw' lc 4 ti "sph"' >> $tg1
    cat $tgr $tg1 >> tmp.gnu
    gnuplot < tmp.gnu
    cat tmp.$mode $tg1 >> $vario_ver

    fln=$fln_vv
    tg1=tmp1
    out='"'$fln'.ps"'
    tit='"Vertical variogram fit ('$mode')"'
    
    echo "set out $out" > tmp.gnu
    echo "set label 1 $tit $labscr" >> tmp.gnu

    echo 'unset y2tics' > $tg1
    echo 'unset y2label' >> $tg1
    echo 'plot \' >> $tg1
    echo '"'$fln'" u 1:3 '$lw' lc 1 ti "smo",\' >> $tg1
    echo '"'$fln'" u 1:4 '$lw' lc 2 ti "exp",\' >> $tg1
    echo '"'$fln'" u 1:5 '$lw' lc 3 ti "gau",\' >> $tg1
    echo '"'$fln'" u 1:6 '$lw' lc 4 ti "sph"' >> $tg1
    cat $tgr $tg1 >> tmp.gnu
    gnuplot < tmp.gnu
    cat tmp.$mode $tg1 >> $vario_hor

#    rm tmp*

    cd $cur
done
# MULTIPLOTS!!
#
font='"Arial,18"' # general font and its size of the title and axes
# set derived gnuplot variables
term="postscript enhanced color font $font"
#term="pngcairo enhanced color crop font $font"
#key="inside bot horizontal font $lfw $lfws"
key="at screen 1,1 Right horizontal font $font $lfws width 0.3 height 0.2"

echo 'set st da l' > $tgr
echo 'set grid' >> $tgr
echo "set term $term" >> $tgr
#
echo 'set xtics nomirror' >> $tgr
echo 'set ytics nomirror' >> $tgr
echo 'set ylab offset 1.5 "(1-x/true) [%]"' >> $tgr
echo "set $psw" >> $tgr
echo "set key $key"  >> $tgr

echo 'set y2tics nomirror' > tmp1
echo 'set y2label offset -1.5 "RMS [%]"' >> tmp1
echo 'set log y2' >> tmp1
echo 'set out "l1_rms_glob.ps"'  >> tmp1
echo 'set origin 0,0' >> tmp1
echo 'set multiplot layout 3,1' >> tmp1
echo 'set log y' >> tmp1


cat $tgr tmp1 $rms_glo > tmp
gnuplot < tmp


echo 'set out "vario_glob.ps"'  > tmp1
echo 'set origin 0,0' >> tmp1
echo 'set multiplot layout 3,1' >> tmp1
cat $tgr tmp1 $vario_glo > tmp
gnuplot < tmp


echo 'set out "vario_hor_ver.ps"'  > tmp1
echo 'set origin 0,0' >> tmp1
echo 'set multiplot layout 3,2' >> tmp1
cat $tgr tmp1 $vario_hor $vario_ver > tmp
gnuplot < tmp

echo 'set out "vario_ver.ps"'  > tmp1
echo 'set origin 0,0' >> tmp1
echo 'set multiplot layout 3,1' >> tmp1
cat $tgr tmp1 $vario_ver > tmp
gnuplot < tmp

exit

#   echo 'set multiplot' >> $tg1
#   echo 'set size 0.52,0.52' >> $tg1 # equal sizes..
#
# setting specific values 
#
# top left, global error
#   echo 'set origin 0,0.5' >> $tg1

    
# variogram plots
echo 'unset y2tics' >> $tg1
echo 'unset y2label' >> $tg1
echo 'unset log y' >> $tg1
echo 'set yrange [*:*]' >> $tg1
echo 'unset key' >> $tg1
# top right global variogram error 
#    echo 'set origin 0.5,0.5' >> tmp.gnu
tit="Global variogram ($mode)"
echo 'set tit "'$tit'"' >> $tg1
echo 'plot \' >> $tg1
echo '"'$fln_vr'" u 1:3 '$lw' lc 1 ti "smo",\' >> $tg1
echo '"'$fln_vr'" u 1:4 '$lw' lc 2 ti "exp",\' >> $tg1
echo '"'$fln_vr'" u 1:5 '$lw' lc 3 ti "gau",\' >> $tg1
echo '"'$fln_vr'" u 1:6 '$lw' lc 4 ti "sph"' >> $tg1

# bottom left, horizontal variogram error
#    echo 'set origin 0,0' >> tmp.gnu
tit="Horizontal variogram ($mode)"
echo 'set tit "'$tit'"' >> $tg1
echo 'plot \' >> $tg1
echo '"'$fln_vh'" u 1:3 '$lw' lc 1 ti "smo",\' >> $tg1
echo '"'$fln_vh'" u 1:4 '$lw' lc 2 ti "exp",\' >> $tg1
echo '"'$fln_vh'" u 1:5 '$lw' lc 3 ti "gau",\' >> $tg1
echo '"'$fln_vh'" u 1:6 '$lw' lc 4 ti "sph"' >> $tg1

# bottom right vertical variogram error
#    echo 'set origin 0.5,0' >> tmp.gnu
tit="Vertical variogram ($mode)"
echo 'set tit "'$tit'"' >> $tg1
echo 'plot \' >> $tg1
echo '"'$fln_vv'" u 1:3 '$lw' lc 1 ti "smo",\' >> $tg1
echo '"'$fln_vv'" u 1:4 '$lw' lc 2 ti "exp",\' >> $tg1
echo '"'$fln_vv'" u 1:5 '$lw' lc 3 ti "gau",\' >> $tg1
echo '"'$fln_vv'" u 1:6 '$lw' lc 4 ti "sph"' >> $tg1
exit    
gnuplot < $tg1
my_pscrop $out
mv `echo $out|sed 's/\"//g'|sed 's/\.ps/\.pdf/g'` ..

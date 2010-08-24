#!/bin/bash

echo running $0 with PID $$ at `uname -n`
date

cur=`pwd`

modes='Exp Gau Sph'


for mode in $modes;do

    if [ -d $mode ];then
	echo working on $mode directory
    else
	echo $mode does not exists.. exiting $0
	exit
    fi

    cd $cur/$mode
    cur2=`pwd`
    fln_glob=$cur2/$mode'_global_l1.dat'
    fln_rms=$cur2/$mode'_rms.dat'
    fln_vr=$cur2/$mode'_variogram_l1_r.dat'
    fln_vh=$cur2/$mode'_variogram_l1_h.dat'
    fln_vv=$cur2/$mode'_variogram_l1_v.dat'
    
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
	    echo sorry, but there are now results to plot
	    echo run $0 without argument
	    exit
	fi
    fi
    
    
# get glob maxvals l1diff (yrange):
    ymaxv=5.e2
    if [ -n $1 ];then
	ymaxv=$1
    fi
    ymax=`awk '!/#/{print $3}' $fln_glob| sort -g|tail -n 1`
    b=`awk '!/#/{print $4}' $fln_glob| sort -g|tail -n 1`
    ymax=`echo $ymax $b|awk '{if($1>$2){print $1}else {print $2}}'`
    b=`awk '!/#/{print $5}' $fln_glob| sort -g|tail -n 1`
    ymax=`echo $ymax $b|awk '{if($1>$2){print $1}else {print $2}}'`
    b=`awk '!/#/{print $6}' $fln_glob| sort -g|tail -n 1`
    ymax=`echo $ymax $b|awk '{if($1>$2){printf("%.0f",$1)}else {printf("%.0f",$2)}}'`
# hier koennte man nochmal abfangen..
    echo $ymax
#ymax=`echo $ymax|awk -v max=$ymaxv '{if($1>max){print max}else{print}}'`
    ymax=`echo $ymax|awk -v max=$ymaxv '{if($1>max){print max}else{print}}'`
    echo $ymax
# get glob minvals l1diff (yrange):
    ymin=`awk '!/#/{print $3}' $fln_glob| sort -g -r|tail -n 1`
    b=`awk '!/#/{print $4}' $fln_glob| sort -g -r|tail -n 1`
    ymin=`echo $ymin $b|awk '{if($1<$2){print $1}else {print $2}}'`
    b=`awk '!/#/{print $5}' $fln_glob| sort -g -r|tail -n 1`
    ymin=`echo $ymin $b|awk '{if($1<$2){print $1}else {print $2}}'`
    b=`awk '!/#/{print $6}' $fln_glob| sort -g -r|tail -n 1`
    ymin=`echo $ymin $b|awk '{if($1<$2){printf("%.0f",$1)}else {printf("%.0f",$2)}}'`
    
# get glob maxvals rms (y2range):
    y2maxv=10
    y2minv=0.8
    y2max=`awk '!/#/{print $3}' $fln_rms| sort -g|tail -n 1`
    b=`awk '!/#/{print $4}' $fln_rms| sort -g|tail -n 1`
    y2max=`echo $y2max $b|awk '{if($1>$2){print $1}else {print $2}}'`
    b=`awk '!/#/{print $5}' $fln_rms| sort -g|tail -n 1`
    y2max=`echo $y2max $b|awk '{if($1>$2){print $1}else {print $2}}'`
    b=`awk '!/#/{print $6}' $fln_rms| sort -g|tail -n 1`
    y2max=`echo $y2max $b|awk '{if($1>$2){printf("%.1f",$1)}else {printf("%.1f",$2)}}'`
# hier koennte man nochmal abfangen..
#y2max=`echo $y2max|awk -v max=$y2maxv '{if($1>max){print max}else{print}}'`
# get glob minvals rms:
    y2min=`awk '!/#/{print $3}' $fln_rms| sort -g -r|tail -n 1`
    b=`awk '!/#/{print $4}' $fln_rms| sort -g -r|tail -n 1`
    y2min=`echo $y2min $b|awk '{if($1<$2){print $1}else {print $2}}'`
    b=`awk '!/#/{print $5}' $fln_rms| sort -g -r|tail -n 1`
    y2min=`echo $y2min $b|awk '{if($1<$2){print $1}else {print $2}}'`
    b=`awk '!/#/{print $6}' $fln_rms| sort -g -r|tail -n 1`
    y2min=`echo $y2min $b|awk '{if($1<$2){printf("%.1f",$1)}else {printf("%.1f",$2)}}'`
# tailoring y2range (min rms..)
    y2min=`echo $y2min|awk '{print $1-0.2}'`
    
    
    echo $ymax $ymin $y2max $y2min
    
    
# declare font size for the legends
    font='"Arial,14"' # general font and its size of the title and axes
# major plots
    lw='w lp lt 1 lw 3.5'   # linetype and width of the major plots
    lfw='"Times,12"' # legend font
    lfw=$font
    lfws='spacing 0.8 samplen 0.4'
    pw='w p pt 1'
    lpw='w l lt 6 lw 2.5'
    psw='pointsize 0.5' # pointsize
    tg1='tmp.gnu'
    tg2=$mode'_l1_evaluations'
    
# set derived gnuplot variables
    term="postscript enhanced color font $font"
#term="pngcairo enhanced color solid transparent crop font $font"
    out='"'$tg2'.ps"'
#key="inside bot horizontal font $lfw $lfws"
    key="below font $lfw $lfws"
    
# setting global gnuplot parameters
    echo 'set st da l' > $tg1
    echo 'set grid' >> $tg1
    echo "set term $term" >> $tg1
# for global model diff plot plot also the RMS
    echo 'set y2tics' >> $tg1
    echo 'set y2tics nomirror' >> $tg1
    echo 'set y2label offset -3 "RMS [%]"' >> $tg1
    echo 'set log y2' >> $tg1
    echo "set y2range [$y2min:$y2max]" >> $tg1
#
    echo 'set xtics nomirror' >> $tg1
    echo 'set xlab offset 0,0.5 "Integral scale /[m]"' >> $tg1
    echo 'set ytics nomirror' >> $tg1
#echo 'set log y' >> $tg1
    echo "set yrange[$ymin:$ymax]" >> $tg1
    echo 'set ylab offset 2 "L_1 difference (1-x/true) [%]"' >> $tg1
    
    echo "set key $key"  >> $tg1
    echo "set $psw" >> $tg1
    echo "set out $out" >> $tg1
    echo 'set multiplot' >> $tg1
    echo 'set size 0.52,0.52' >> $tg1 # equal sizes..
#
# setting specific values 
#
# top left, global error
    echo 'set origin 0,0.5' >> $tg1
    tit="Global model ($mode)"
    echo 'set tit "'$tit'"' >> $tg1
    echo 'plot \' >> $tg1
    echo '"'$fln_glob'" u 1:3 '$lw' lc 1 ti "smo",\' >> $tg1
    echo '"'$fln_rms'" u 1:3 axes x1y2 '$lpw' lc 1 ti "(rms)",\' >> $tg1
    echo '"'$fln_glob'" u 1:4 '$lw' lc 2 ti "exp",\' >> $tg1
    echo '"'$fln_rms'" u 1:4 axes x1y2 '$lpw' lc 2 ti "(rms)",\' >> $tg1
    echo '"'$fln_glob'" u 1:5 '$lw' lc 3 ti "gau",\' >> $tg1
    echo '"'$fln_rms'" u 1:5 axes x1y2 '$lpw' lc 3 ti "(rms)",\' >> $tg1
    echo '"'$fln_glob'" u 1:6 '$lw' lc 4 ti "sph",\' >> $tg1
    echo '"'$fln_rms'" u 1:6 axes x1y2 '$lpw' lc 4 ti "(rms)"' >> $tg1
    
# variogram plots
    echo 'unset y2tics' >> $tg1
    echo 'unset y2label' >> $tg1
    echo 'unset log y' >> $tg1
    echo 'set yrange [*:*]' >> $tg1
    echo 'unset key' >> $tg1
# top right global variogram error 
    echo 'set origin 0.5,0.5' >> tmp.gnu
    tit="Global variogram ($mode)"
    echo 'set tit "'$tit'"' >> $tg1
    echo 'plot \' >> $tg1
    echo '"'$fln_vr'" u 1:3 '$lw' lc 1 ti "smo",\' >> $tg1
    echo '"'$fln_vr'" u 1:4 '$lw' lc 2 ti "exp",\' >> $tg1
    echo '"'$fln_vr'" u 1:5 '$lw' lc 3 ti "gau",\' >> $tg1
    echo '"'$fln_vr'" u 1:6 '$lw' lc 4 ti "sph"' >> $tg1
    
# bottom left, horizontal variogram error
    echo 'set origin 0,0' >> tmp.gnu
    tit="Horizontal variogram ($mode)"
    echo 'set tit "'$tit'"' >> $tg1
    echo 'plot \' >> $tg1
    echo '"'$fln_vh'" u 1:3 '$lw' lc 1 ti "smo",\' >> $tg1
    echo '"'$fln_vh'" u 1:4 '$lw' lc 2 ti "exp",\' >> $tg1
    echo '"'$fln_vh'" u 1:5 '$lw' lc 3 ti "gau",\' >> $tg1
    echo '"'$fln_vh'" u 1:6 '$lw' lc 4 ti "sph"' >> $tg1
    
# bottom right vertical variogram error
    echo 'set origin 0.5,0' >> tmp.gnu
    tit="Vertical variogram ($mode)"
    echo 'set tit "'$tit'"' >> $tg1
    echo 'plot \' >> $tg1
    echo '"'$fln_vv'" u 1:3 '$lw' lc 1 ti "smo",\' >> $tg1
    echo '"'$fln_vv'" u 1:4 '$lw' lc 2 ti "exp",\' >> $tg1
    echo '"'$fln_vv'" u 1:5 '$lw' lc 3 ti "gau",\' >> $tg1
    echo '"'$fln_vv'" u 1:6 '$lw' lc 4 ti "sph"' >> $tg1
    
    gnuplot < $tg1
    my_pscrop $out
    mv `echo $out|sed 's/\"//g'|sed 's/\.ps/\.pdf/g'` ..
    cd $cur
done
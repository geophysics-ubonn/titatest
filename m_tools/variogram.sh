#!/bin/bash

# declare font size for the legends
font='"Helvetica" 16' # general font and its size of the title and axes
# major plots
lw='w l lw 4'   # linetype and width of the major plots
lfw='Times,12' # legend font
lfws='spacing 1.1'
psw='pointsize 1.5' # pointsize

fac='100.0' # factor for the rel diff

if [ -z $3 ];then
    variofunc='Experimental/Variogram function,'
else
    if [ "$3" == "Exp" ];then
	variofunc='{/Symbol g}(h)=va*(1-EXP(-(3h/a))),'
    elif [ "$3" == "Gau" ];then
	variofunc='{/Symbol g}(h)=va*(1-EXP(-(3h/a)^2)),'
    elif [ "$3" == "Sph" ];then
	variofunc='{/Symbol g}(h)=va((1.5(h/a)-.5(h/a)**3),1),'
    else
	variofunc='Experimental/Variogram function,'
    fi
fi

if [ -z $2 ];then
    echo trying to deduce arguments
    cur=`pwd`
    bla=`basename $cur`
    ax=`echo $bla|tr '_' ' '|awk '{print $1}'`
    ay=`echo $bla|tr '_' ' '|awk '{print $2}'`
    if [[ -z $ax || -z $ay ]];then
	echo "usage $0 with arguments ax and ah";
	exit
    fi
else
    ax=$1
    ay=$2
fi

ah=`echo $ax $ay | awk '{print $1/$2}' `

echo "$0 $variofunc with integral scale horizontal:$ax, ratio $ah"

tg1='tmp.gnu'
tg2='varioplots.tex'

echo 'set st da l' > $tg1
echo 'set grid' >> $tg1
echo "set term pos enh col sol $font" >> $tg1
echo 'set xtics nomirror' >> $tg1
echo 'set ytics nomirror' >> $tg1
echo 'set y2tics' >> $tg1
echo 'set y2tics nomirror' >> $tg1
echo 'set log y2' >> $tg1
echo 'set xlab offset 0,0.5 "Lag h [m]"' >> $tg1
echo 'set ylab offset 2,0 "sv(h)=1/2N(h) {/Symbol S}_i(Z(m_i+h)-Z(m_i))^2"' >> $tg1
echo 'set y2la offset -1.5 "Rel. diff [%]"' >> $tg1
echo 'se key inside right bot nobox noreverse font "'$lfw'"' $lfws >> $tg1
echo "set $psw" >> $tg1
tit="ax=$ax, ax/az=$ah"
echo $tit
# r-variogram
echo 'set tit "'$variofunc' '$tit', h = (x^2+z^2)^{1/2}"' >> $tg1
echo 'set out "variograms.ps"' >> $tg1

true=`tail -n 1 true/inv.variogram | awk '{print $3}' `
smo=`tail -n 1 smo/inv.variogram | awk '{print $3}' `
exp=`tail -n 1 exp/inv.variogram | awk '{print $3}' `
gau=`tail -n 1 gau/inv.variogram | awk '{print $3}' `
sph=`tail -n 1 sph/inv.variogram | awk '{print $3}' `

echo $ax $ay $true $smo $exp $gau $sph > variogram_l1_r.dat

echo "true=$true" >> $tg1
echo "smo=$smo" >> $tg1
echo "exp=$exp" >> $tg1
echo "gau=$gau" >> $tg1
echo "sph=$sph" >> $tg1

paste true/inv.variogram smo/inv.variogram|awk -v t=$true -v b=$smo '!/#/{a=(1-($6/$2)*(t/b));printf("%f\t%f\t%f\t%f\n",$1,sqrt(a*a)*100,$2,$6)}' > smo/diffs_inv.dat
paste true/inv.variogram exp/inv.variogram|awk -v t=$true -v b=$exp '!/#/{a=(1-($6/$2)*(t/b));printf("%f\t%f\t%f\t%f\n",$1,sqrt(a*a)*100,$2,$6)}' > exp/diffs_inv.dat
paste true/inv.variogram gau/inv.variogram|awk -v t=$true -v b=$gau '!/#/{a=(1-($6/$2)*(t/b));printf("%f\t%f\t%f\t%f\n",$1,sqrt(a*a)*100,$2,$6)}' > gau/diffs_inv.dat
paste true/inv.variogram sph/inv.variogram|awk -v t=$true -v b=$sph '!/#/{a=(1-($6/$2)*(t/b));printf("%f\t%f\t%f\t%f\n",$1,sqrt(a*a)*100,$2,$6)}' > sph/diffs_inv.dat

l1_smo=`awk '{sum+=$2} END {printf("%.2f\n",sum/(NR))}' smo/diffs_inv.dat`
l1_exp=`awk '{sum+=$2} END {printf("%.2f\n",sum/(NR))}' exp/diffs_inv.dat`
l1_gau=`awk '{sum+=$2} END {printf("%.2f\n",sum/(NR))}' gau/diffs_inv.dat`
l1_sph=`awk '{sum+=$2} END {printf("%.2f\n",sum/(NR))}' sph/diffs_inv.dat`

echo 'plot \' >> $tg1
echo '"true/inv.variogram" u 1:($3/true) '$lw' lc 0 ti "{/Symbol g}(h)",\' >> $tg1
echo '"true/inv.variogram" u 1:($2/true) '$lw' lc 1 ti "sv(h),true",\' >> $tg1
echo '"exp/inv.variogram" u 1:($2/exp) '$lw' lc 2 ti "sv(h),exp",\' >> $tg1
echo '"gau/inv.variogram" u 1:($2/gau) '$lw'  lc 3 ti "sv(h),gau",\' >> $tg1
echo '"sph/inv.variogram" u 1:($2/sph) '$lw'  lc 4 ti "sv(h),sph",\' >> $tg1
echo '"smo/inv.variogram" u 1:($2/smo) '$lw'  lc 5 ti "sv(h),smo",\' >> $tg1
# factor true/x is from the scaling: $exp/$true*(true/exp) -> normalized 
echo '"exp/diffs_inv.dat" u 1:($2*(true/exp)) axes x1y2 w p lc 2 ti "exp/true(L1='$l1_exp' %)",\' >> $tg1
echo '"gau/diffs_inv.dat" u 1:($2*(true/gau)) axes x1y2 w p lc 3 ti "gau/true(L1='$l1_gau' %)",\' >> $tg1
echo '"sph/diffs_inv.dat" u 1:($2*(true/sph)) axes x1y2 w p lc 4 ti "sph/true(L1='$l1_sph' %)",\' >> $tg1
echo '"smo/diffs_inv.dat" u 1:($2*(true/smo)) axes x1y2 w p lc 5 ti "smo/true(L1='$l1_smo' %)"' >> $tg1


# x-variogram
tit="$ax"
echo 'set tit "'$variofunc' a = ax ='$tit', h = hx = x"' >> $tg1
echo 'set out "variograms_x.ps"' >> $tg1

true=`tail -n 1 true/inv.variogram_x | awk '{print $3}' `
smo=`tail -n 1 smo/inv.variogram_x | awk '{print $3}' `
exp=`tail -n 1 exp/inv.variogram_x | awk '{print $3}' `
gau=`tail -n 1 gau/inv.variogram_x | awk '{print $3}' `
sph=`tail -n 1 sph/inv.variogram_x | awk '{print $3}' `

echo $ax $ay $true $smo $exp $gau $sph > variogram_l1_h.dat

echo "true=$true" >> $tg1
echo "smo=$smo" >> $tg1
echo "exp=$exp" >> $tg1
echo "gau=$gau" >> $tg1
echo "sph=$sph" >> $tg1


paste true/inv.variogram_x smo/inv.variogram_x|awk -v t=$true -v b=$smo '!/#/{a=(1-($6/$2)*(t/b));printf("%f\t%f\t%f\t%f\n",$1,sqrt(a*a)*100,$2,$6)}' > smo/diffs_inv_x.dat
paste true/inv.variogram_x exp/inv.variogram_x|awk -v t=$true -v b=$exp '!/#/{a=(1-($6/$2)*(t/b));printf("%f\t%f\t%f\t%f\n",$1,sqrt(a*a)*100,$2,$6)}' > exp/diffs_inv_x.dat
paste true/inv.variogram_x gau/inv.variogram_x|awk -v t=$true -v b=$gau '!/#/{a=(1-($6/$2)*(t/b));printf("%f\t%f\t%f\t%f\n",$1,sqrt(a*a)*100,$2,$6)}' > gau/diffs_inv_x.dat
paste true/inv.variogram_x sph/inv.variogram_x|awk -v t=$true -v b=$sph '!/#/{a=(1-($6/$2)*(t/b));printf("%f\t%f\t%f\t%f\n",$1,sqrt(a*a)*100,$2,$6)}' > sph/diffs_inv_x.dat

l1_smo=`awk '{sum+=$2} END {printf("%.2f\n",sum/(NR))}' smo/diffs_inv_x.dat`
l1_exp=`awk '{sum+=$2} END {printf("%.2f\n",sum/(NR))}' exp/diffs_inv_x.dat`
l1_gau=`awk '{sum+=$2} END {printf("%.2f\n",sum/(NR))}' gau/diffs_inv_x.dat`
l1_sph=`awk '{sum+=$2} END {printf("%.2f\n",sum/(NR))}' sph/diffs_inv_x.dat`

echo 'plot \' >> $tg1
echo '"true/inv.variogram_x" u 1:($3/true) '$lw' lc 0 ti "{/Symbol g}(h)",\' >> $tg1
echo '"true/inv.variogram_x" u 1:($2/true) '$lw' lc 1 ti "sv(h),true",\' >> $tg1
echo '"exp/inv.variogram_x" u 1:($2/exp) '$lw' lc 2 ti "sv(h),exp",\' >> $tg1
echo '"gau/inv.variogram_x" u 1:($2/gau) '$lw'  lc 3 ti "sv(h),gau",\' >> $tg1
echo '"sph/inv.variogram_x" u 1:($2/sph) '$lw'  lc 4 ti "sv(h),sph",\' >> $tg1
echo '"smo/inv.variogram_x" u 1:($2/smo) '$lw'  lc 5 ti "sv(h),smo",\' >> $tg1
echo '"exp/diffs_inv_x.dat" u 1:($2) axes x1y2 w p lc 2 ti "exp/true(L1='$l1_exp' %)",\' >> $tg1
echo '"gau/diffs_inv_x.dat" u 1:($2) axes x1y2 w p lc 3 ti "gau/true(L1='$l1_gau' %)",\' >> $tg1
echo '"sph/diffs_inv_x.dat" u 1:($2) axes x1y2 w p lc 4 ti "sph/true(L1='$l1_sph' %)",\' >> $tg1
echo '"smo/diffs_inv_x.dat" u 1:($2) axes x1y2 w p lc 5 ti "smo/true(L1='$l1_smo' %)"' >> $tg1

# y-variogram
tit=`echo $ax $ah | awk '{print $1/$2}' `
echo 'set tit "'$variofunc' a = az ='$tit', h = hz = z"' >> $tg1
echo 'set out "variograms_y.ps"' >> $tg1

true=`tail -n 1 true/inv.variogram_y  | awk '{print $3}' `
smo=`tail -n 1 smo/inv.variogram_y  | awk '{print $3}' `
exp=`tail -n 1 exp/inv.variogram_y  | awk '{print $3}' `
gau=`tail -n 1 gau/inv.variogram_y  | awk '{print $3}' `
sph=`tail -n 1 sph/inv.variogram_y  | awk '{print $3}' `

echo $ax $ay $true $smo $exp $gau $sph > variogram_l1_v.dat

echo "true=$true" >> $tg1
echo "smo=$smo" >> $tg1
echo "exp=$exp" >> $tg1
echo "gau=$gau" >> $tg1
echo "sph=$sph" >> $tg1


paste true/inv.variogram_y smo/inv.variogram_y|awk -v t=$true -v b=$smo '!/#/{a=(1-($6/$2)*(t/b));printf("%f\t%f\t%f\t%f\n",$1,sqrt(a*a)*100,$2,$6)}' > smo/diffs_inv_y.dat
paste true/inv.variogram_y exp/inv.variogram_y|awk -v t=$true -v b=$exp '!/#/{a=(1-($6/$2)*(t/b));printf("%f\t%f\t%f\t%f\n",$1,sqrt(a*a)*100,$2,$6)}' > exp/diffs_inv_y.dat
paste true/inv.variogram_y gau/inv.variogram_y|awk -v t=$true -v b=$gau '!/#/{a=(1-($6/$2)*(t/b));printf("%f\t%f\t%f\t%f\n",$1,sqrt(a*a)*100,$2,$6)}' > gau/diffs_inv_y.dat
paste true/inv.variogram_y sph/inv.variogram_y|awk -v t=$true -v b=$sph '!/#/{a=(1-($6/$2)*(t/b));printf("%f\t%f\t%f\t%f\n",$1,sqrt(a*a)*100,$2,$6)}' > sph/diffs_inv_y.dat
l1_smo=`awk '{sum+=$2} END {printf("%.2f\n",sum/(NR))}' smo/diffs_inv_y.dat`
l1_exp=`awk '{sum+=$2} END {printf("%.2f\n",sum/(NR))}' exp/diffs_inv_y.dat`
l1_gau=`awk '{sum+=$2} END {printf("%.2f\n",sum/(NR))}' gau/diffs_inv_y.dat`
l1_sph=`awk '{sum+=$2} END {printf("%.2f\n",sum/(NR))}' sph/diffs_inv_y.dat`


echo 'plot \' >> $tg1
echo '"true/inv.variogram_y" u 1:($3/true) '$lw' lc 0 ti "{/Symbol g}(h)",\' >> $tg1
echo '"true/inv.variogram_y" u 1:($2/true) '$lw' lc 1 ti "sv(h),true",\' >> $tg1
echo '"exp/inv.variogram_y" u 1:($2/exp) '$lw' lc 2 ti "sv(h),exp",\' >> $tg1
echo '"gau/inv.variogram_y" u 1:($2/gau) '$lw'  lc 3 ti "sv(h),gau",\' >> $tg1
echo '"sph/inv.variogram_y" u 1:($2/sph) '$lw'  lc 4 ti "sv(h),sph",\' >> $tg1
echo '"smo/inv.variogram_y" u 1:($2/smo) '$lw'  lc 5 ti "sv(h),smo",\' >> $tg1
echo '"exp/diffs_inv_y.dat" u 1:($2) axes x1y2 w p lc 2 ti "exp/true(L1='$l1_exp' %)",\' >> $tg1
echo '"gau/diffs_inv_y.dat" u 1:($2) axes x1y2 w p lc 3 ti "gau/true(L1='$l1_gau' %)",\' >> $tg1
echo '"sph/diffs_inv_y.dat" u 1:($2) axes x1y2 w p lc 4 ti "sph/true(L1='$l1_sph' %)",\' >> $tg1
echo '"smo/diffs_inv_y.dat" u 1:($2) axes x1y2 w p lc 5 ti "smo/true(L1='$l1_smo' %)"' >> $tg1

gnuplot < $tg1

my_pscrop ./

echo '\documentclass[12pt,a4]{article}' > $tg2
echo '\usepackage{amsmath,amssymb,mathrsfs}' >> $tg2
echo '\usepackage[T1]{fontenc}' >> $tg2
echo '\usepackage[utf8]{inputenc}' >> $tg2
echo '\usepackage{listings}' >> $tg2
echo '\lstset{basicstyle=\tiny}' >> $tg2
echo '\usepackage{pstricks,epsf,graphicx,psfrag,epsfig}' >> $tg2
echo '\usepackage{epstopdf,boxedminipage,fancybox,subfigure,float}' >> $tg2
echo '\usepackage{fancyhdr}' >> $tg2
echo '\oddsidemargin-15mm' >> $tg2
echo '\parindent0mm' >> $tg2
echo '\parskip0mm' >> $tg2
echo '\textheight24cm' >> $tg2
echo '\textwidth19cm' >> $tg2
echo '\unitlength1mm' >> $tg2
echo '\topmargin-2cm' >> $tg2
echo '\headheight0cm' >> $tg2
echo '\footskip3cm' >> $tg2
echo '\begin{document}' >> $tg2
echo '\pagestyle{fancy}' >> $tg2
echo '\fancyhead[C]{\Large\bf CRTomo inversion results}' >> $tg2
echo '\includegraphics[width=.5\textwidth]{true/model.pdf}' >> $tg2
echo '\includegraphics[width=.5\textwidth]{smo/model.pdf}' >> $tg2
echo '\includegraphics[width=.5\textwidth]{exp/model.pdf}' >> $tg2
echo '\includegraphics[width=.5\textwidth]{gau/model.pdf}' >> $tg2
echo '\includegraphics[width=.5\textwidth]{sph/model.pdf}' >> $tg2
echo '\clearpage' >> $tg2
if [ -e diffs_exp_modl.pdf ];then
    echo 'found diff plots!!'
    echo '\fancyhead[C]{\Large\bf Relative differences}' >> $tg2    
    for x in diffs_*_modl.pdf;do
	echo including $x
	echo "\includegraphics[width=.5\textwidth]{$x}" >> $tg2
    done
    echo '\clearpage' >> $tg2
fi
echo '\fancyhead[C]{\Large\bf Corresponding Variograms}' >> $tg2
echo '\includegraphics[width=.5\textwidth]{variograms.pdf}' >> $tg2
echo '\includegraphics[width=.5\textwidth]{variograms_x.pdf}' >> $tg2
echo '\includegraphics[width=.5\textwidth]{variograms_y.pdf}' >> $tg2
echo '\end{document}' >> $tg2


pdflatex $tg2

echo 'finished'

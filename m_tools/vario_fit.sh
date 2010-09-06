#!/bin/bash

gplt=tmp.gnu

if [ -n "$1" ];then
    echo copying $1 to variogram.dat
else

    if [ -e variogram.dat ];then
	echo found variogram.dat
    else
	echo no input nor variogram.dat found
	exit
    fi
fi


echo 've(x)=sexp * (1 - exp( -(3*x/aexp) ) )' > $gplt
echo 'fit ve(x) "variogram.dat" u 1:2 via sexp,aexp' >> $gplt
# make LSQ-fit
gnuplot < $gplt >& variofit_ve.dat
# get the output
sexp=`grep -A 4 "Final set" variofit_ve.dat|awk '/sexp/{printf("%.3f\n",sqrt($3*$3))}'`
dsexp=`grep -A 4 "Final set" variofit_ve.dat|awk '/sexp/{printf("%.3f\n",$5)}'`
aexp=`grep -A 4 "Final set" variofit_ve.dat|awk '/aexp/{printf("%.2f\n",sqrt($3*$3))}'`
daexp=`grep -A 4 "Final set" variofit_ve.dat|awk '/aexp/{printf("%.2f\n",$5)}'`
# prepare plot title
texp="c (1 - exp( -(3h/a) )); c=$sexp{/Symbol \261}$dsexp, a=$aexp{/Symbol \261}$daexp\n"

echo 'vg(x)=sgau * (1 - exp( -(3*x/agau)**2.) )' > $gplt
echo 'fit vg(x) "variogram.dat" u 1:2 via sgau,agau' >> $gplt
gnuplot < $gplt >& variofit_vg.dat
sgau=`grep -A 4 "Final set" variofit_vg.dat|awk '/sgau/{printf("%.3f\n",sqrt($3*$3))}'`
dsgau=`grep -A 4 "Final set" variofit_vg.dat|awk '/sgau/{printf("%.3f\n",$5)}'`
agau=`grep -A 4 "Final set" variofit_vg.dat|awk '/agau/{printf("%.2f\n",sqrt($3*$3))}'`
dagau=`grep -A 4 "Final set" variofit_vg.dat|awk '/agau/{printf("%.2f\n",$5)}'`
# prepare plot title
tgau="c (1 - exp( -(3h/a)^2 )); c=$sgau{/Symbol \261}$dsgau, a=$agau{/Symbol \261}$dagau"

echo 'vs(x)=ssph * (1.5 * (x/asph) - .5 * (x/asph)**3 )' > $gplt
echo 'fit vs(x) "variogram.dat" u 1:2 via ssph,asph' >> $gplt
gnuplot < $gplt >& variofit_vs.dat
ssph=`grep -A 4 "Final set" variofit_vs.dat|awk '/ssph/{printf("%.3f\n",sqrt($3*$3))}'`
dssph=`grep -A 4 "Final set" variofit_vs.dat|awk '/ssph/{printf("%.3f\n",$5)}'`
asph=`grep -A 4 "Final set" variofit_vs.dat|awk '/asph/{printf("%.2f\n",sqrt($3*$3))}'`
dasph=`grep -A 4 "Final set" variofit_vs.dat|awk '/asph/{printf("%.2f\n",$5)}'`

tsph="c (1.5(h/a) - .5(h/a)^3); c=$ssph{/Symbol \261}$dssph, a=$asph{/Symbol \261}$dasph\n"

echo '# Exp Gau Sph / Errors (absolute)' > variofit.dat
echo $aexp $agau $asph $daexp $dagau $dasph >> variofit.dat

gplt=$gplt'2'


echo 'set tit "Variogram fit"' > $gplt
echo 'set xla offset 0,.5 "Lag (h)"' >> $gplt
echo 'set yla offset 1.5 "Variogram"' >> $gplt
echo 'set grid' >> $gplt
echo 'set term pos enh col sol 20' >> $gplt
echo 'set pointsize 1.4' >> $gplt
echo 'set out "variogram_fit.ps"' >> $gplt
echo 'set key bot right Right font "Arial,18" samplen .3 spacing 2.5' >> $gplt
echo "ve(x)=$sexp*(1 - exp(-(3*x/$aexp)))" >> $gplt
echo "vg(x)=$sgau*(1 - exp(-(3*x/$agau)**2))" >> $gplt
echo "vs(x)=$ssph*(1.5*(x/$asph) - .5*(x/$asph)**3)" >> $gplt
echo 'p\'>> $gplt
echo '"variogram.dat" w p lc 0 pt 7 ti "Semivariogram",\' >> $gplt
echo 've(x) w l lw 6 lc 1 ti "'$texp'",\' >> $gplt
echo 'vg(x) w l lw 6 lc 2 ti "'$tgau'",\' >> $gplt
echo 'vs(x) w l lw 6 lc 3 ti "'$tsph'"' >> $gplt

gnuplot < $gplt 

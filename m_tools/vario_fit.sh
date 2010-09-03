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
echo 'vg(x)=sgau * (1 - exp( -(3*x/agau)**2.) )' >> $gplt
echo 'vs(x)=ssph * (1.5 * (x/asph) - .5 * (x/asph)**3 )' >> $gplt
echo 'fit ve(x) "variogram.dat" u 1:2 via sexp,aexp' >> $gplt
echo 'fit vg(x) "variogram.dat" u 1:2 via sgau,agau' >> $gplt
echo 'fit vs(x) "variogram.dat" u 1:2 via ssph,asph' >> $gplt
echo 'save var "variofit.dat"' >> $gplt

gnuplot < $gplt

gplt=$gplt'2'
sgau=`awk '/sgau =/{printf("%.3f\n",$3)}' variofit.dat`
agau=`awk '/agau =/{printf("%.3f\n",$3)}' variofit.dat`
sexp=`awk '/sexp =/{printf("%.3f\n",$3)}' variofit.dat`
aexp=`awk '/aexp =/{printf("%.3f\n",$3)}' variofit.dat`
ssph=`awk '/ssph =/{printf("%.3f\n",$3)}' variofit.dat`
asph=`awk '/asph =/{printf("%.3f\n",$3)}' variofit.dat`

tgau="$sgau(1-exp(-(3h/$agau)**2))"
texp="$sexp(1-exp(-(3h/$aexp)))"
tsph="$ssph(1.5(h/$asph)-.5(h/$asph)**3)"

echo 'set tit "Variogram least squares fit"' > $gplt
echo 'set xla offset 0,.5 "Lag (h)"' >> $gplt
echo 'set yla offset 1 "Variogram"' >> $gplt
echo 'set term pos enh col sol 20' >> $gplt
echo 'set pointsize 1.5' >> $gplt
echo 'set out "variogram_fit.eps"' >> $gplt
echo 'set key bot right' >> $gplt
echo "ve(x)=$sexp*(1 - exp(-(3*x/$aexp)))" >> $gplt
echo "vg(x)=$sgau*(1 - exp(-(3*x/$agau)**2))" >> $gplt
echo "vs(x)=$ssph*(1.5*(x/$asph) - .5*(x/$asph)**3)" >> $gplt
echo 'p\'>> $gplt
echo '"variogram.dat" w p lc 0,\' >> $gplt
echo 'vg(x) w l lw 3 lc 1 ti "'$tgau'",\' >> $gplt
echo 've(x) w l lw 3 lc 2 ti "'$texp'",\' >> $gplt
echo 'vs(x) w l lw 3 lc 3 ti "'$tsph'"' >> $gplt

gnuplot < $gplt 

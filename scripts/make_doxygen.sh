#!/bin/bash
commit=`git log | head -1|awk '{print $2}'`
branch=`git branch|awk '/\*/{print $2}'`
mydate=`date +%c|sed 's/\ /-/g'`

if [ -e doxy.inp ];then
	echo doxy.inp exists.
else
	doxygen -g doxy.inp
fi
tag="$1"
awk -v tag=$tag '{if($1=="PROJECT_NAME"){print $1" = "tag}else{print}}' doxy.inp > tmp.inp
tag=$branch
awk -v tag=$tag '{if($1=="PROJECT_NUMBER"){print $1" = "tag}else{print}}' tmp.inp > doxy.inp
outdir="$2"
awk -v tag=$outdir '{if($1=="OUTPUT_DIRECTORY"){print $1" = "tag}else{print}}' tmp.inp > doxy.inp

my_list_of_values="OPTIMIZE_FOR_FORTRAN EXTRACT_ALL GENERATE_TREE_VIEW HAVE_DOT CALL_GRAPH CALLER_GRAPH GENERATE_LATEX USE_PDFLATEX PDF_HYPERLINKS LATEX_SOURCE_CODE GENERATE_DOCBOOK GENERATE_TREEVIEW"
for member in $my_list_of_values;do

	awk -v member=$member '{if($1==member){print $1" = Yes"}else{print}}' doxy.inp > tmp.inp
	mv tmp.inp doxy.inp

done
doxygen doxy.inp
cd $outdir/latex
make

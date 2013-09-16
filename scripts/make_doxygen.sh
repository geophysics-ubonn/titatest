#!/bin/bash
commit=$(git log | head -1|awk '{print $2}')
branch=$(git branch|awk '/\*/{print $2}')
mydate=$(date +%c|sed 's/\ /-/g')

if [ -f dox.inp ];then
	echo "taking existent doxy.inp"
else
	echo "creating default doxy input file"
	doxygen -g doxy.inp
fi

tag1="$1"
tag2="$branch"
tag3="$2"
echo "edit PROJECT_NAME = $tag1"
echo "edit PROJECT_NUMBER = $tag2"
echo "edit OUTPUT_DIR = $tag3"
awk -v t1="$tag1" -v t2="$tag2" -v t3="$tag3" '{ 
						if($1 == "PROJECT_NAME" ){print $1" = "t1}
						else if($1 == "PROJECT_NUMBER" ){print $1" = "t2}
						else if($1 == "OUTPUT_DIR" ){print $1" = "t3}
						else{print}
						}' doxy.inp > tmp.inp

cp tmp.inp dox.inp
my_list_of_values="OPTIMIZE_FOR_FORTRAN EXTRACT_ALL GENERATE_TREE_VIEW HAVE_DOT CALL_GRAPH CALLER_GRAPH GENERATE_LATEX USE_PDFLATEX PDF_HYPERLINKS LATEX_SOURCE_CODE GENERATE_DOCBOOK GENERATE_TREEVIEW"
#my_list_of_values="OPTIMIZE_FOR_FORTRAN EXTRACT_ALL GENERATE_TREE_VIEW HAVE_DOT CALL_GRAPH CALLER_GRAPH GENERATE_LATEX USE_PDFLATEX PDF_HYPERLINKS LATEX_SOURCE_CODE GENERATE_DOCBOOK GENERATE_TREEVIEW"
for member in $my_list_of_values;do

	awk -v member=$member '{if($1==member){print $1" = Yes"}else{print}}' doxy.inp > tmp.inp
	mv tmp.inp doxy.inp

done
doxygen doxy.inp
cd $tag3/latex
make

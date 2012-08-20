#!/bin/bash
cur_branch=`git branch | awk '/\*/{print $2}'`
echo "current branch: $cur_branch"
liste=`git branch|sed 's/\*//g' `
echo $liste
for x in $liste;do
	git checkout $x
	git pull
done
git checkout $cur_branch

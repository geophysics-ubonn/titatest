#!/bin/bash
echo "CHARACTER (256) :: my_git_version(3)" > my_git_version.h
commit=`git log | head -1|awk '{print $2}'`
branch=`git branch|awk '/\*/{print $2}'`
mydate=`date +%c|sed 's/\ /-/g'`
echo 'my_git_version(1) = "'$branch'"' >> my_git_version.h
echo 'my_git_version(2) = "'$commit'"' >> my_git_version.h
echo 'my_git_version(3) = "'$mydate'"' >> my_git_version.h

#!/bin/bash
# store git status (branch, commit, compiler, os, date) into my_git_version.h.
# This file will be compiled into the binaries.
# Run in src/

# get the character specification
#grep "CHARACTER (256)" get_git_ver.f90 | sed 's/,PUBLIC//g' | sed 's/version/my_git_version/g' > my_git_version.h
commit=$( git log | head -1 | awk '{print $2}' )
branch=$( git branch | awk '/\*/{print $2}' )
mydate=$( date +%c | sed 's/\ /-/g' )
compiler=$1
myos=$( uname -o )
echo 'Character(256) :: my_git_version(5)' > my_git_version.h
echo 'my_git_version(1) = "'$branch'"' >> my_git_version.h
echo 'my_git_version(2) = "'$commit'"' >> my_git_version.h
echo 'my_git_version(3) = "'$mydate'"' >> my_git_version.h
echo 'my_git_version(4) = "'$compiler'"' >> my_git_version.h
echo 'my_git_version(5) = "'$myos'"' >> my_git_version.h


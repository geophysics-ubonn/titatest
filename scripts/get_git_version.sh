#!/bin/bash
# get the character specification
grep "CHARACTER (256)" get_git_ver.f90 | sed 's/,PUBLIC//g' | sed 's/version/my_git_version/g' > my_git_version.h
commit=$( git log | head -1 | awk '{print $2}' )
branch=$( git branch | awk '/\*/{print $2}' )
mydate=$( date +%c | sed 's/\ /-/g' )
compiler=$1
myos=$( uname -o )
echo 'my_git_version(1) = "'$branch'"' >> my_git_version.h
echo 'my_git_version(2) = "'$commit'"' >> my_git_version.h
echo 'my_git_version(3) = "'$mydate'"' >> my_git_version.h
echo 'my_git_version(4) = "'$compiler'"' >> my_git_version.h
echo 'my_git_version(5) = "'$myos'"' >> my_git_version.h


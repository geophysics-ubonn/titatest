#!/bin/bash
# change version number in those places:
#	* configure.ac
#	* debian/changelog

# the Debian-changelog can be changed with:
# debchange --newversion 1.6 --distribution stable

# generate configure, configure
# the Debian package always has the same prefix
# export FCFLAGS="-O4 -march=core2 -ftree-vectorize -ffast-math -funroll-loops -finline-functions -fopenmp"
CRT_PREFIX="dev" ./autogen.sh

# create archive
make dist-bzip2

# delete if necessary
test -d package && rm -r package

mkdir package

mv crtomomod*.bz2 package
cd package
tar xvjf crtomomod*.bz2
cd crtomomod*
dpkg-buildpackage -us -uc
# echo "Now enter package/crtomomod*/ and type dpkg-buildpackage"

#!/bin/sh
# Run this script in a fresh git cloned directory to generate the configure
# script.

# generate src/Makefile.am
cd src
../scripts/make_crerr.sh
cd ..

# generate cutmckee/Makefile.am
cd cutmckee
../scripts/make_crerr.sh
cd ..

# generate git credentials
cd src
../scripts/get_git_version.sh
cd ..

# create AUTHORS file
echo "AK - A. Kemna (kemna@geo.uni-bonn.de)" > AUTHORS
echo "JK - J. Kenkel (jkenkel@geo.uni-bonn.de)" >> AUTHORS
echo "RM - R. Martin (rmartin@geo.uni-bonn.de)" >> AUTHORS
echo "MW - M. Weigand (mweigand@geo.uni-bonn.de)" >> AUTHORS

touch NEWS ChangeLog

echo -n "aclocal... "
aclocal

# run autoreconf multiple times to detect missing files
echo -n "autoreconf (1)... "
autoreconf -v
# add missing files
echo -n "automake... "
automake --add-missing
automake

# run again, now that we added all files
echo -n "autoreconf (2)... "
autoreconf -v

./configure

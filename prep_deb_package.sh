#!/bin/bash

# generate configure
./prep_autotools.sh

# configure and create archive
./configure
make dist-bzip2

# delete if necessary
if [ -d package ]; then
	rm -r package
fi

mkdir package

mv crtomomod*.bz2 package
cd package
tar xvjf crtomomod*.bz2


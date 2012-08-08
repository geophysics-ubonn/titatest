#!/bin/bash

# generate configure, configure
./autogen.sh

# create archive
make dist-bzip2

# delete if necessary
if [ -d package ]; then
	rm -r package
fi

mkdir package

mv crtomomod*.bz2 package
cd package
tar xvjf crtomomod*.bz2


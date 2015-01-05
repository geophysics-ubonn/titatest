#!/bin/sh
make clean

rm curmckee/Makefile.am
rm man/Makefile.am
rm minimalbeispiele//Makefile.am
rm src/crerror.h
rm cutmckee/crerror.h
rm -r config

# Die Makefiles werden beim Wechseln der branches überschrieben und gelöscht, also müssen eventuelle Objekte gelöscht werden
rm src/*.o
rm src/*.mod
rm cutmckee/*.o
rm cutmckee/*.mod

# crerror.h wird automatisch generiert
rm inc/crerror.h
rm cutmck/crerror.h

# Lösche alle von autotools generierten Dateien

# Autotools
rm Makefile src/Makefile cutmckee/Makefile
rm Makefile.in src/Makefile.in cutmckee/Makefile.in
rm man/Makefile

rm aclocal.m4

rm config.guess
rm config.log
rm config.status
rm config.sub

rm Man/Makefile

rm configure

rm INSTALL
rm install-sh
rm libtool
rm ltmain.sh
rm missing
test -z depcomp || rm depcomp

rm -r autom4te.cache

rm -rf package

rm -rf crtomomod-1.0/
rm -f debian/crtomomod.debhelper.log

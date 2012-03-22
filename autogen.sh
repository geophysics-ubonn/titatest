#!/bin/sh

# Diese Dateien müssen wir noch schreiben. Bis dahin werden sie immer leer erstellt (und gelöscht)
echo "AK - A. Kemna (kemna@geo.uni-bonn.de)" > AUTHORS
echo "JK - J. Kenkel (jkenkel@geo.uni-bonn.de)" >> AUTHORS
echo "RM - R. Martin (rmartin@geo.uni-bonn.de)" >> AUTHORS
echo "MW - M. Weigand (mweigand@geo.uni-bonn.de)" >> AUTHORS

CRT_cfg=man/crtomo.cfg

if [ -e$CRT_cfg ];then
	cp $CRT_cfg README
else
	touch README
fi	

touch NEWS ChangeLog

aclocal
# Das erste autoreconf bemerkt, welche Dateien fehlen
autoreconf
# Wenn möglich, werden Standarddateien kopiert
automake --add-missing
# Jetzt sollten alle Dateien vorhanden sein, und das configure-Skript kann erstellt werden
autoreconf

#!/bin/sh

# Diese Dateien müssen wir noch schreiben. Bis dahin werden sie immer leer erstellt (und gelöscht)
echo "AK - A. Kemna (kemna@geo.uni-bonn.de)" > AUTHORS
echo "JK - J. Kenkel (jkenkel@geo.uni-bonn.de)" >> AUTHORS
echo "RM - R. Martin (rmartin@geo.uni-bonn.de)" >> AUTHORS
echo "MW - M. Weigand (mweigand@geo.uni-bonn.de)" >> AUTHORS

touch NEWS ChangeLog

echo -n "aclocal... "
aclocal
# Das erste autoreconf bemerkt, welche Dateien fehlen
echo -n "autoreconf (1)... "
autoreconf
# Wenn möglich, werden Standarddateien kopiert
echo -n "automake... "
automake --add-missing
# Jetzt sollten alle Dateien vorhanden sein, und das configure-Skript kann erstellt werden
echo -n "autoreconf (2)... "
autoreconf


exit
# now invoke configure
echo "configure"
./configure
# Wenn möglich, werden Standarddateien kopiert
echo -n "automake... "
automake --add-missing
echo "ok"

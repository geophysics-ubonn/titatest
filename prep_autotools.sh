#!/bin/sh

# Diese Dateien müssen wir noch schreiben. Bis dahin werden sie immer leer erstellt (und gelöscht)
touch NEWS README AUTHORS ChangeLog

# Das erste autoreconf bemerkt, welche Dateien fehlen
autoreconf
# Wenn möglich, werden Standarddateien kopiert
automake --add-missing
# Jetzt sollten alle Dateien vorhanden sein, und das configure-Skript kann erstellt werden
autoreconf

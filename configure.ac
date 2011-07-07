
AC_PREREQ([2.64])
# das Projekt initialisieren, versionsnummer und emailadresse angeben
AC_INIT([crtomomod], [1.0], [mweigand@geo.uni-bonn.de])
# eine beliebige Sourcedatei angeben, damit der Sourceordner bekannt ist
AC_CONFIG_SRCDIR([src/inv.f90])


# Automake soll aufgerufen werden (nicht nur autoconf)
AM_INIT_AUTOMAKE

# CRTomo and CRMod error generation
AC_CONFIG_COMMANDS_PRE([cd src; ./../make_crerr.sh; cd ..;])

# CutMcK error generation
AC_CONFIG_COMMANDS_PRE([cd cutmckee; ./../make_crerr.sh; cd ..;])

# old (development) Makefile should be moved to a save place
#AC_CONFIG_COMMANDS_PRE([mv src/Makefile src/old.Makefile])
#AC_CONFIG_COMMANDS_PRE([mv CutMcK/Makefile CutMcK/old.Makefile])

# den autotools mitteilen, das Fortran benutzt wird. Mit folgenden Aufrufen werden die
# entsprechenden Variablen initialisiert 
AC_PROG_FC
AC_PROG_F77

# Lapack und Blas vorhanden ? Checks für configure
#ACX_LAPACK

# An diesen Stellen sollen Makefiles von configure erstellt werden
AC_CONFIG_FILES([Makefile
		 src/Makefile
		 cutmckee/Makefile
		 man/Makefile])

AC_OUTPUT

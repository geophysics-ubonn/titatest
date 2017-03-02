AC_PREREQ(2.63)
# das Projekt initialisieren, versionsnummer und emailadresse angeben
AC_INIT([CRTomoMod], [1.6], [mweigand@geo.uni-bonn.de])
# eine beliebige Sourcedatei angeben, damit der Sourceordner bekannt ist
AC_CONFIG_SRCDIR([src/inv.f90])
AC_CONFIG_AUX_DIR([config])

mySCRIPTS=scripts

# apparently the call to AC_PROG_FF/AC_PROG_F77 will change the user-supplied
# settings for the compilers. Therefore, save and restore the content later

# normally we should not set FCFLAGS and FFLAGS (these are for user-supplied
# settings). However, we want to change some machine specific stuff for the
# Debian package

default_flags=" -O4 -march=native -ftree-vectorize -ffast-math -funroll-loops -finline-functions "
test -z "${FCFLAGS}" && FCFLAGS="${default_flags}"
test -z "${FFLAGS}" && FFLAGS="${default_flags}"

FCFLAGS_SAVE="${FCFLAGS}"
FFLAGS_SAVE="${FFLAGS}"

# default compiler flags
# AM_FC_FLAGS=" -fopenmp "
AM_FC_FLAGS=" "
AC_SUBST(AM_FC_FLAGS)

# den autotools mitteilen, das Fortran benutzt wird. Mit folgenden Aufrufen werden die
# entsprechenden Variablen initialisiert
AC_PROG_FC([gfortran])
AC_PROG_F77([gfortran])

# restore user compiler flags
FCFLAGS="${FCFLAGS_SAVE}"
FFLAGS="${FFLAGS_SAVE}"

AC_PREFIX_DEFAULT($HOME)

# Automake soll aufgerufen werden (nicht nur autoconf)
AM_INIT_AUTOMAKE

# CRTomo and CRMod error generation
AC_CONFIG_COMMANDS_PRE([cd src; ./../$mySCRIPTS/make_crerr.sh ; cd ..;])

# CRTomo git version generation
AC_CONFIG_COMMANDS_PRE([cd src; ./../$mySCRIPTS/get_git_version.sh $myCompiler ; cd ..;])

# CutMcK error generation
AC_CONFIG_COMMANDS_PRE([cd cutmckee; ./../$mySCRIPTS/make_crerr.sh ; cd ..;])


# old (development) Makefile should be moved to a save place
#AC_CONFIG_COMMANDS_PRE([mv src/Makefile src/old.Makefile])
#AC_CONFIG_COMMANDS_PRE([mv CutMcK/Makefile CutMcK/old.Makefile])


# Lapack und Blas vorhanden ? Checks für configure
#ACX_LAPACK

# An diesen Stellen sollen Makefiles von configure erstellt werden
AC_CONFIG_FILES([Makefile
		 src/Makefile
		 cutmckee/Makefile
		 man/Makefile])
# generiere output
AC_OUTPUT

# Minimalbeispiel precision arithemtics
.SUFFIX:	.o .f90

P1		= minimal_prec
P2		= minimal_omp


bin_PROGRAMS	= $(P1) $(P2)

FFLAG90         = -O4 -march=native -ftree-vectorize -ffast-math -funroll-loops -finline-functions

FFLAG90         = -pg -fbounds-check -Wuninitialized -O -ftrapv -pedantic \
		-fdump-fortran-original  -fdump-fortran-optimized \
		-ffpe-trap=list  -fdump-core -fbacktrace -fdump-parse-tree \
		-fimplicit-none -fno-signed-zeros -ffinite-math-only -fopenmp

LDFLAGS		= -lm # link math lib

P1_SRC		= $(P1).f90
P2_SRC		= $(P2).f90

minimal_prec_SOURCES	=	$(P1_SRC)
minimal_prec_FCFLAGS	=	$(FFLAG90)
minimal_prec_LDFLAGS	=	$(LDFLAGS)

minimal_omp_SOURCES	=	$(P2_SRC)
minimal_omp_FCFLAGS	=	$(FFLAG90)
minimal_omp_LDFLAGS	=	$(LDFLAGS)


clean-local:	
		if [ `ls -1  | sed '/\.mod/!d' | wc -l` -gt 0 ]; then  rm *.mod; fi

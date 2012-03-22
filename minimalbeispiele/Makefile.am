# Minimalbeispiel precision arithemtics
.SUFFIX:	.o .f90

P1		= minimal_prec
P2		= minimal_omp


bin_PROGRAMS	= $(P1) $(P2)

AM_FC_FLAGS = -O4 -march=native -ftree-vectorize -ffast-math -funroll-loops -finline-functions

AM_DBG_FLAGS  = -pg -fbounds-check -Wuninitialized -O -ftrapv -pedantic \
		-fdump-fortran-original  -fdump-fortran-optimized \
		-ffpe-trap=list  -fdump-core -fbacktrace -fdump-parse-tree \
		-fimplicit-none -fno-signed-zeros -ffinite-math-only -fopenmp

AM_LD_FLAGS = -lm # link math lib

P1_SRC		= tic_toc.f90 minimal-prec.f90
P2_SRC		= minimal-omp.f90

minimal_prec_SOURCES	=	$(P1_SRC)
minimal_prec_FCFLAGS	=	$(AM_FC_FLAGS)
minimal_prec_LDFLAGS	=	$(AM_LD_FLAGS)

minimal_omp_SOURCES	=	$(P2_SRC)
minimal_omp_FCFLAGS	=	$(AM_FC_FLAGS) -fopenmp
minimal_omp_LDFLAGS	=	$(AM_LD_FLAGS)


clean-local:	
		if [ `ls -1  | sed '/\.mod/!d' | wc -l` -gt 0 ]; then  rm *.mod; fi

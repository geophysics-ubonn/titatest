# $Id: Makefile,v 1.0 2009/07/02 12:43:39 roland Exp $
#
#      Makefile for CutMck (default)
#       

RM		= rm -f
CP		= cp -f
MV		= mv -f
WPATH 		= ~/bin

F90		= gfortran
F77		= gfortran
FFLAG90         = -O3 -march=native -ftree-vectorize -fexpensive-optimizations -ffast-math
#FFLAG90         = -C -Wall -mcmodel=medium
FFLAGMPI        = -I/usr/include/lam
FFLAGMPI        = 
FLIBMPI         = -L/usr/lib/lam/lib -llammpio -llamf77mpi -lmpi -llam -lutil -ldl -lnsl
FLIBMPI         = 
FLIBLINBLAS     = -llapack -lblas
FLIBLINBLAS     = 
FLIB            = -lm
FLIBF77         = -lm

# das hier chek obs ein bin im home gibt
C1		= cbn
# name des Programms 
PROG	        = CutMcK
# target
PR1		= ctm

# default targets
all:		$(PR1)

.SILENT:	cbn
#################################################
#
OBJ		= cutmck.o relem.o welem.o relectr.o welectr.o \
		  get_error.o
# suffixes
.SUFFIXES:	.f90 .f .c .o .mod
#.SUFFIXES:	.f90 .f .c .o
# rules
.f90.mod:		
		$(F90) $(FFLAG90) $(FFLAGMPI) -c $<

.mod.o:		
		$(F90) $(FFLAG90) $(FFLAGMPI) -c $<

.f90.o:		
		$(F90) $(FFLAG90) $(FFLAGMPI) -c $<

.f.o:		
		$(F90) $(FFLAG90) $(FFLAGMPI) -c $<

.c.o:		
		$(CC) $(CFLAGS) $(CFLAGMPI) -c $<
# dependencies	
get_error.o:	error.txt get_error.f90
		./make_crerr.sh
		$(F90) $(FFLAG90) -c get_error.f90

# targets

cbn:		
		echo "Pruefe ~/bin"
		if [ -d ~/bin ]; then \
			echo "ok"; \
		else \
			echo "Du hast kein bin in deinem home.--"; \
			mkdir ~/bin; \
		fi

ctm:		$(OBJ) 
		$(F90) $(FFLAG90) $(FFLAGMPI) -o $(PROG) \
		$(OBJ) \
		$(FLIBLINBLAS) $(FLIBMPI) $(FLIB)

clean:		
		$(RM) *.o *~ *.vo V* mod.* *.mod *.flc *.d \
			$(PROG)

install:	$(PROG)
		$(CP) $(PROG) $(WPATH)

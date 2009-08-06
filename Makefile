# $Id: Makefile,v 1.0 2009/07/02 12:43:39 roland Exp $
#
#      Makefile for CRTomo (default)
#       

RM		= rm -f
CP		= cp -f
MV		= mv -f
WPATH 		= ~/bin

F90		= gfortran
F77		= gfortran
FFLAG90         = -O3 -march=native -ftree-vectorize -fexpensive-optimizations -ffast-math
FFLAG90         = -C
FFLAGMPI        = -I/usr/include/lam
FFLAGMPI        = 
FLIBMPI         = -L/usr/lib/lam/lib -llammpio -llamf77mpi -lmpi -llam -lutil -ldl -lnsl
FLIBMPI         = 
FLIBLINBLAS     = -llapack -lblas
#FLIBLINBLAS     = 
FLIB            = -lm
FLIBF77         = -lm

PR1		= crt

PR2		= crm

all:		$(PR1) $(PR2)

#.SILENT:	all crt crm
# default targets

crt:		*.for inv.f
		$(F90) $(FFLAG90) $(FFLAGMPI) $(FLIBLINBLAS) -o CRTomo \
		*.for inv.f
		$(CP) CRTomo $(WPATH)

crm:		*.for fem.f
		$(F90) $(FFLAG90) $(FFLAGMPI) $(FLIBLINBLAS) -o CRMod \
		*.for fem.f
		$(CP) CRMod $(WPATH)
install:		
		$(CP) CRTomo $(WPATH)
		$(CP) CRMod $(WPATH)

clean:		
		$(RM) CRTomo CRmod *~ *.mod

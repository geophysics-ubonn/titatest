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
FLIBLINBLAS     = 
FLIB            = -lm
FLIBF77         = -lm

PR1		= crt

PR2		= crm

linux:		all

# default targets

crt:		*.for
		if [ -e "fem.for" ];then $(MV) fem.for cr.fem; fi
		$(F90) $(FFLAG90) $(FFLAGMPI) -o CRTomo *.for
		$(CP) CRTomo $(WPATH)
		if [ -e "cr.fem" ];then $(MV) cr.fem fem.for; fi

crm:		*.for
		if [ -e "inv.for" ];then $(MV) inv.for cr.inv; fi
		$(F90) $(FFLAG90) $(FFLAGMPI) -o CRMod *.for
		$(CP) CRMod $(WPATH)
		if [ -e "cr.inv" ];then $(MV) cr.inv inv.for; fi


all:		
		crt crm

install:		
		$(CP) CRTomo $(WPATH)
		$(CP) CRMod $(WPATH)

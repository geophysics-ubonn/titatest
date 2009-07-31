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

CI		= checkinv 

IC 		= invcheck

CF		= checkfem

FC 		= femcheck

all:		$(PR1) $(PR2)

#.SILENT:	all crt crm
# default targets

crt:		*.for
		if [ -e "fem.for" ];then \
			echo $(MV) fem.for cr.fem;\
			$(MV) fem.for cr.fem;\
	        fi
		$(F90) $(FFLAG90) $(FFLAGMPI) -o CRTomo *.for
		$(CP) CRTomo $(WPATH)
		if [ -e "cr.fem" ];then \
			echo $(MV) cr.fem fem.for; \
			$(MV) cr.fem fem.for; \
		fi

crm:		*.for
		if [ -e "inv.for" ];then \
			echo $(MV) inv.for cr.inv; \
			$(MV) inv.for cr.inv; \
		fi
		$(F90) $(FFLAG90) $(FFLAGMPI) -o CRMod *.for
		$(CP) CRMod $(WPATH)
		if [ -e "cr.inv" ];then \
			echo $(MV) cr.inv inv.for; \
			$(MV) cr.inv inv.for; \
		fi
install:		
		$(CP) CRTomo $(WPATH)
		$(CP) CRMod $(WPATH)

clean:		
		$(RM) CRTomo CRmod

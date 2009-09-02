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

PR3		= mtools

f90crt		= get_unit.o make_noise.o alloci.o get_error.o

f90crm		= alloci.o

all:		$(PR1) $(PR2) $(PR3)


# rules
#.f90.mod:		
#		$(F90) $(FFLAG90) -c $<
#
#.mod.o:		
#		$(F90) $(FFLAG90) -c $<
#

#.f90.o:		
#		$(F90) $(FFLAG90) -c $<

#.SILENT:	all crt crm
# default targets
################################## F90 targets..
make_noise.o:	make_noise.f90
		$(F90) $(FFLAG90) -c make_noise.f90

get_error.o:	get_error.f90
		$(F90) $(FFLAG90) -c get_error.f90

get_unit.o:	get_unit.f90
		$(F90) $(FFLAG90) -c get_unit.f90

alloci.o:	alloci.f90
	        $(F90) $(FFLAG90) -c alloci.f90
###################################

crt:		*.for inv.f $(f90crt)
		$(F90) $(FFLAG90) $(FFLAGMPI) $(FLIBLINBLAS) -o CRTomo \
		*.for inv.f $(f90crt)
		$(CP) CRTomo $(WPATH)

crm:		*.for fem.f $(f90crt)
		$(F90) $(FFLAG90) $(FFLAGMPI) $(FLIBLINBLAS) -o CRMod \
		*.for fem.f $(f90crt)
		$(CP) CRMod $(WPATH)

mtools:		
		$(CP) m_tools/crtomo_plot.sh $(WPATH)
		$(CP) m_tools/crtomo_run.sh $(WPATH)
		$(CP) m_tools/plot_cur_crmod $(WPATH)
		$(CP) m_tools/plotCRTmod_batch.m $(WPATH)

install:		
		$(CP) CRTomo $(WPATH)
		$(CP) CRMod $(WPATH)

clean:		
		$(RM) CRTomo CRMod *~ *.mod *.o m_tools/*~

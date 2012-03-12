# $Id: Makefile,v 1.0 2009/07/02 12:43:39 roland Exp $
#
#      Makefile for CutMck (default)
#       

RM		= rm -f
CP		= cp -f
MV		= mv -f

F90		= gfortran
F77		= gfortran
FFLAG90         = -O4 -march=native -ftree-vectorize -fexpensive-optimizations -ffast-math
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
# machine name
BRANCH		= $(shell git branch|awk '/\*/{print $$2}')
MACHINE		= $(shell uname -n)
PREFIX		= $(BRANCH)

# query for OS dependant flags
ifeq    ($(OS), Msys)
	WPATH           = "C:\\MinGW\\bin"
else
	WPATH           = $${HOME}/bin
	PREFIX          = $(BRANCH)_$(MACHINE)_$(F90)
endif
					
# default targets
all:		$(PR1) install

.SILENT:	cbn
#################################################
#
OBJ		= cutmck.o relem.o welem.o relectr.o welectr.o

OBJF90		= elemmod.o electrmod.o

OBJERR		= get_error.o
###################################################
$(OBJ):		%.o : %.f	
		$(F90) $(FFLAG90) -c $<

$(OBJF90):	%.o : %.f90	
		$(F90) $(FFLAG90) -c $<
# dependencies	
$(OBJERR):	error.txt get_error.f90
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

ctm:		$(OBJF90) $(OBJ) $(OBJERR)
		$(F90) $(FFLAG90) $(FFLAGMPI) -o $(PROG)_$(PREFIX) \
		$(OBJ) $(OBJF90) $(OBJERR) \
		$(FLIBLINBLAS) $(FLIBMPI) $(FLIB)

clean:		
		$(RM) *.o *~ *.vo V* mod.* *.mod *.flc *.d \
		$(PROG)_$(PREFIX)

install:	$(PROG)_$(PREFIX)
		$(CP) $(PROG)_$(PREFIX) $(WPATH)

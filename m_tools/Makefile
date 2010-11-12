# $Id: Makefile,v 1.0 2009/07/02 12:43:39 roland Exp $
#
#      Makefile for CRTomo (default)
#       

RM		= rm -f
CP		= cp -f
MV		= mv -f
WPATH 		= ~/bin
###############################
C1		= cbn
# kopiert die matlab tools
PRM		= mtools

all:		$(PRM)

.SILENT:	cbn
# default targets

SCRIPTS	=	crtomo_*.sh crmod_*.sh plt_*_crt* plot_*.sh plot_cur_* \
		plotCRTmod_batch*.m vario*.sh stochastic_*.sh

cbn:		
		echo "Pruefe ~/bin"
		if [ -d ~/bin ]; then \
			echo "ok"; \
		else \
			echo "Du hast kein bin in deinem home.--"; \
			mkdir ~/bin; \
		fi

mtools:		$(C1)
		$(CP) $(SCRIPTS) $(WPATH)

install:	$(C1)
		$(CP) $(SCRIPTS) $(WPATH)

clean:		
		$(RM) *~

chmo:
		echo "chmod " $(SCRIPTS)
		chmod 0755 $(SCRIPTS)

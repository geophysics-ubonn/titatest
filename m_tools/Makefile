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

SCRIPTS	=	crtomo_plot.sh crtomo_run.sh plt_all_crt plot_cur_crtomo \
		plot_cur_crmod plot_cur_crt plt_inv_ctr plotCRTmod_batch.m \
		variogram.sh vario_fit.sh

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

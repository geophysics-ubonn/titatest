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
################################## F90 targets..
cbn:		
		echo "Pruefe ~/bin"
		if [ -d ~/bin ]; then \
			echo "ok"; \
		else \
			echo "Du hast kein bin in deinem home.--"; \
			mkdir ~/bin; \
		fi

mtools:		$(C1)
		$(CP) crtomo_plot.sh $(WPATH)
		$(CP) crtomo_run.sh $(WPATH)
		$(CP) plt_all_crt $(WPATH)
		$(CP) plot_cur_crmod $(WPATH)
		$(CP) plot_cur_crtomo $(WPATH)
		$(CP) plot_cur_crt $(WPATH)
		$(CP) plotCRTmod_batch.m $(WPATH)

install:	$(C1)
		$(CP) crtomo_plot.sh $(WPATH)
		$(CP) crtomo_run.sh $(WPATH)
		$(CP) plt_all_crt $(WPATH)
		$(CP) plot_cur_crmod $(WPATH)
		$(CP) plot_cur_crtomo $(WPATH)
		$(CP) plot_cur_crt $(WPATH)
		$(CP) plotCRTmod_batch.m $(WPATH)

clean:		
		$(RM) *~

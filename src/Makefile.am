
## determine arch type and add 64 bit directives, if possible
#ARCH           =$(shell uname -a)


###
bin_PROGRAMS = CRTomo CRMod
###
### CRMod source files
f90crm		= alloci.f90 femmod.f90 datmod.f90 \
		  invmod.f90 sigmamod.f90 electrmod.f90 modelmod.f90 \
		  elemmod.f90 wavenmod.f90 randbmod.f90 errmod.f90 konvmod.f90 \
		  pathmod.f90 ompmod.f90

fcrm		= fem.f90

f90crmsub	= bbsens.f90 bessi0.f90 bessi1.f90 bessk0.f90 bessk1.f90 \
		  beta.f90 bkfak.f90 bpot.f90 \
		  bsendc.f90 bsens.f90 bsensi.f90 \
		  bsytop.f90 bvolt.f90 bvolti.f90 chareal.f90 \
		  chkpol.f90 choldc.f90 chol.f90 \
		  elem1.f90 elem3.f90 elem4.f90 \
		  elem5.f90 elem8.f90 filpat.f90 \
		  kompadc.f90 kompbdc.f90 kompb.f90 \
		  potana.f90 precal.f90 \
		  randb2.f90 randb.f90 randdc.f90 rdatm.f90 \
		  relectr.f90 relem.f90 rrandb.f90 rsigma.f90 refsig.f90 \
		  rtrafo.f90 rwaven.f90 scalab.f90 scaldc.f90 sort.f90 \
		  vredc.f90 vre.f90 wdatm.f90 \
		  wkpot.f90 wout.f90 wpot.f90 wsens.f90 \
		  gammln.f90 gaulag.f90 gauleg.f90 intcha.f90 kompab.f90 \
		  tic_toc.f90 make_noise.f90 get_unit.f90 get_git_ver.f90 get_error.f90 
		  


# CRTomo objects
f90crt		= alloci.f90 femmod.f90 \
		  datmod.f90 invmod.f90 modelmod.f90 cjgmod.f90 sigmamod.f90 \
		  electrmod.f90  elemmod.f90 wavenmod.f90 \
		  randbmod.f90 errmod.f90 konvmod.f90 pathmod.f90 \
		  invhpmod.f90 ompmod.f90 cg_mod.f90 tic_toc.f90 variomodel.f90 bsmatm_mod.f90 bmcm_mod.f90 \
		  brough_mod.f90 get_git_ver.f90

f90crtsub	= bbsedc.f90 bbsens.f90 besp_elem.f90 \
		  bessi0.f90 bessi1.f90 bessk0.f90 bessk1.f90 \
		  beta.f90 bkfak.f90 blam0.f90 bnachbar.f90 bpar.f90 bpot.f90 \
		  bsendc.f90 bsens.f90 bsensi.f90 \
		  bsigm0.f90 bsigma.f90 bsytop.f90 bvolt.f90 bvolti.f90 chareal.f90 \
		  chkpol.f90 choldc.f90 chol.f90 dmisft.f90 \
		  elem1.f90 elem3.f90 elem4.f90 \
		  elem5.f90 elem8.f90 filpat.f90 \
		  chold.f90 cholz.f90 linvd.f90 linvz.f90 \
		  mdian1.f90 parfit.f90 potana.f90 precal.f90 make_noise.f90 rall.f90 \
		  gammln.f90 gaulag.f90 gauleg.f90 intcha.f90 kompab.f90 \
		  kompadc.f90 kompbdc.f90 kompb.f90 kont1.f90 kont2.f90 \
		  randb2.f90 randb.f90 randdc.f90 rdati.f90 rdatm.f90 \
		  refsig.f90 relectr.f90 relem.f90 rrandb.f90 rsigma.f90 \
		  rtrafo.f90 rwaven.f90 scalab.f90 scaldc.f90 sort.f90 \
		  update.f90 vredc.f90 vre.f90 wdatm.f90 \
		  wkpot.f90 wout.f90 wpot.f90 wsens.f90 \
		  gauss_dble.f90 gauss_cmplx.f90 get_unit.f90 \
		   bvariogram.f90 get_error.f90
		  

fcrt		= inv.f90


CRTomo_SOURCES = $(f90crt) $(f90crtsub) $(fcrt)
CRTomo_FFLAGS = -O3 -march=native -ftree-vectorize -fexpensive-optimizations -ffast-math -I ../inc 
# check if compiling on 64 bit, add memory model flag
#if "$(ARCH)" == "x86_64"
#	CRTomo_FFLAGS = $(CRTomo_FFLAGS) -mcmodel=medium
#endif


CRTomo_FCFLAGS = $(CRTomo_FFLAGS)
CRTomo_LDFLAGS = -lm

CRMod_SOURCES = $(f90crm) $(fcrm) $(f90crmsub)
CRMod_FFLAGS = $(CRTomo_FFLAGS)
CRMod_FCFLAGS = $(CRMod_FFLAGS)
CRMod_LDFLAGS = $(CRTomo_LDFLAGS)


# die Includedatei crerr.h muss neu generiert werden, wenn error.txt sich ändert
#get_error.f90:	../error.txt
#		cd ..; ./../make_crerr.sh; cd src;


# zusätzliche clean-rule, um Fortranmodule zu löschen, falls sie existieren
clean-local:	
		if [ `ls -1  | sed '/\.mod/!d' | wc -l` -gt 0 ]; then  rm *.mod; fi

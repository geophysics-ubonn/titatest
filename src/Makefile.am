
## determine arch type and add 64 bit directives, if possible
ARCH           =$(shell uname -a)


###
bin_PROGRAMS = CRTomo CRMod
###
### CRMod source files
f90crm		= alloci.o femmod.o datmod.o \
		  invmod.o sigmamod.o electrmod.o modelmod.o \
		  elemmod.o wavenmod.o randbmod.o errmod.o konvmod.o \
		  pathmod.o ompmod.o

fcrm		= fem.o

f90crmsub	= bbsens.o bessi0.o bessi1.o bessk0.o bessk1.o \
		  beta.o bkfak.o bpot.o \
		  bsendc.o bsens.o bsensi.o \
		  bsytop.o bvolt.o bvolti.o chareal.o \
		  chkpol.o choldc.o chol.o \
		  elem1.o elem3.o elem4.o \
		  elem5.o elem8.o filpat.o \
		  kompadc.o kompbdc.o kompb.o \
		  potana.o precal.o \
		  randb2.o randb.o randdc.o rdatm.o \
		  relectr.o relem.o rrandb.o rsigma.o refsig.o \
		  rtrafo.o rwaven.o scalab.o scaldc.o sort.o \
		  vredc.o vre.o wdatm.o \
		  wkpot.o wout.o wpot.o wsens.o \
		  gammln.o gaulag.o gauleg.o intcha.o kompab.o \
		  tic_toc.o make_noise.o get_unit.o


# CRTomo objects
f90crt		= alloci.o femmod.o \
		  datmod.o invmod.o cjgmod.o sigmamod.o \
		  electrmod.o modelmod.o elemmod.o wavenmod.o \
		  randbmod.o errmod.o konvmod.o pathmod.o \
		  invhpmod.o ompmod.o

f90crtsub	= bbsedc.o bbsens.o besp_elem.o \
		  bessi0.o bessi1.o bessk0.o bessk1.o \
		  beta.o bkfak.o blam0.o bnachbar.o bpar.o bpot.o \
		  bsendc.o bsens.o bsensi.o \
		  bsigm0.o bsigma.o bsytop.o bvolt.o bvolti.o chareal.o \
		  chkpol.o choldc.o chol.o dmisft.o \
		  elem1.o elem3.o elem4.o \
		  elem5.o elem8.o filpat.o \
		  chold.o cholz.o linvd.o linvz.o \
		  mdian1.o parfit.o potana.o precal.o rall.o \
		  gammln.o gaulag.o gauleg.o intcha.o kompab.o \
		  kompadc.o kompbdc.o kompb.o kont1.o kont2.o \
		  randb2.o randb.o randdc.o rdati.o rdatm.o \
		  refsig.o relectr.o relem.o rrandb.o rsigma.o \
		  rtrafo.o rwaven.o scalab.o scaldc.o sort.o \
		  update.o vredc.o vre.o wdatm.o \
		  wkpot.o wout.o wpot.o wsens.o \
		  gauss_dble.o gauss_cmplx.o get_unit.o \
		  make_noise.o tic_toc.o variomodel.o bvariogram.o \
		  cg_mod.o bsmatm_mod.o bmcm_mod.o brough_mod.o \

fcrt		= inv.o


CRTomo_SOURCES = $(f90crt) $(fcrt) $(f90crtsub)
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
get_error.f90:	../inc/error.txt
		cd ../inc; ./../make_crerr.sh; cd ../src;


# zusätzliche clean-rule, um Fortranmodule zu löschen, falls sie existieren
clean-local:	
		if [ `ls -1  | sed '/\.mod/!d' | wc -l` -gt 0 ]; then  rm *.mod; fi

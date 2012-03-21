# $Id: Makefile,v 1.0 2009/07/02 12:43:39 roland Exp $
#
#      Makefile for CRTomo (default)
#       

FFLAG90         = -O4 -march=native -ftree-vectorize -ffast-math -funroll-loops -finline-functions
## Only un-comment for debug purpose
#FFLAG90         = -pg -fbounds-check -Wuninitialized -O -ftrapv -pedantic \
		-fdump-fortran-original  -fdump-fortran-optimized \
		-fdump-core -fbacktrace -fdump-parse-tree \
		-fimplicit-none -fno-signed-zeros -ffinite-math-only -fopenmp
############################################################################
# 			ifort compiler flags 				   #
# make sure to comment all other compiler flags from above in order to use # 
# ifort									   #
############################################################################
#F90		= ifort
#FFLAG90		= -O3 -fast -parallel
## Only un-comment for debug purpose
#FFLAG90         = -C -g -debug all -check all -implicitnone \
		-warn unused -fp-stack-check -heap-arrays -ftrapuv \
		-check pointers -check bounds -openmp

bin_PROGRAMS =	CRTomo_testing \
		CRMod_testing

main_mods = 	alloci.f90 \
		datmod.f90 \
		electrmod.f90 elemmod.f90 errmod.f90 \
		femmod.f90 \
		invmod.f90 \
		konvmod.f90 \
		modelmod.f90 \
		ompmod.f90 \
		pathmod.f90 \
		randbmod.f90 \
		sigmamod.f90 \
		wavenmod.f90


main_subs = 	bbsens.f90 bessi0.f90 beta.f90 bessi1.f90 bessk0.f90 \
		bessk1.f90 bkfak.f90 bsytop.f90 \
		chareal.f90 chkpol.f90 chol.f90 \
		elem1.f90 elem3.f90 elem4.f90 elem5.f90 elem8.f90 \
		filpat.f90 \
		gammln.f90 gaulag.f90 gauleg.f90 get_err.f90 \
		intcha.f90 \
		kompab.f90 kompb.f90 \
		make_noise.f90 \
		potana.f90 precal.f90 \
		randb.f90 randb2.f90 rdatm.f90 relectr.f90 relem.f90 \
		rrandb.f90 rsigma.f90 refsig.f90 rtrafo.f90 rwaven.f90 \
		scalab.f90 sort.f90 \
		tic_toc.f90 \
		vre.f90 \
		wdatm.f90 wkpot.f90 wout.f90

# CRTomo main program
crt_main = 	inv.f90

# CRTomo additional modules
crt_mods =	cjgmod.f90 \
		invhpmod.f90 

# CRTomo additional subroutines 
crt_subs = 	bmcm_mod.f90 brough_mod.f90 bsendc.f90 bsensi.f90 bsmatm_mod.f90 \
		bvariogram.f90 bvolti.f90 \
		cg_mod.f90 chold.f90 choldc.f90 cholz.f90 \
		dmisft.f90 \
		gauss_dble.f90 gauss_cmplx.f90 \
		kont1.f90 kont2.f90 kompadc.f90 kompbdc.f90 \
		linvd.f90 linvz.f90 \
		mdian1.f90 \
		parfit.f90 \
		rdati.f90 \
		scaldc.f90 \
		update.f90 \
		variomodel.f90 vredc.f90 \
		wout_up.f90

# CRMod main program
crm_main =	fem.f90
# CRMomo additional subroutines 
crt_subs = 	bpot.f90 bsens.f90 bvolt.f90 \
		wdatm.f90 wpot.f90 wsens.f90


CRTomo_testing_SOURCES = $(main_mods) $(main_subs) $(crt_main) $(crt_subs)
CRTomo_testing_FFLAGS =  $(FFLAGS)
CRTomo_testing_FCFLAGS = $(CRTomo_FFLAGS)
CRTomo_testing_LDFLAGS = -lm

CRMod_testing_SOURCES = $(f90crm) $(fcrm) $(f90crmsub)
CRMod_testing_FFLAGS = $(CRTomo_FFLAGS)
CRMod_testing_FCFLAGS = $(CRMod_FFLAGS)
CRMod_testing_LDFLAGS = $(CRTomo_LDFLAGS)


# clean-rule to erase fortran modules
clean-local:	
		if [ `ls -1  | sed '/\.mod/!d' | wc -l` -gt 0 ]; then  rm *.mod; fi

EXTRA_DIST = 	error.txt \
		my_git_version.h


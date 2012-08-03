bin_PROGRAMS = CutMcK_dev_kreide

cOBJ		= relem.f welem.f relectr.f welectr.f
cOBJF90		= elemmod.f90 electrmod.f90
cOBJERR		= get_error.f90

CutMcK_dev_kreide_SOURCES = $(cOBJF90) $(cOBJ) $(cOBJERR) cutmck.f 
CutMcK_dev_kreide_FFLAGS =  -O3 -march=native -ftree-vectorize -fexpensive-optimizations -ffast-math
CutMcK_dev_kreide_FCFLAGS = $(CutMcK_FFLAGS)
CutMcK_dev_kreide_LDFLAGS = -lm # -llapack -lblas

# zusätzliche clean-rule, um Fortranmodule zu löschen, falls sie existieren
clean-local:
		if [ `ls -1  | sed '/\.mod/!d' | wc -l` -gt 0 ]; then  rm *.mod; fi
		rm -f Makefile.am

EXTRA_DIST = error.txt err.fin # file with error text messages

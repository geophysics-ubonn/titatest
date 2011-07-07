EXTRA_DIST = error.txt

bin_PROGRAMS = CutMcK

cOBJ		= cutmck.o relem.o welem.o relectr.o welectr.o
cOBJF90		= elemmod.o electrmod.o
cOBJERR		= get_error.o

CutMcK_SOURCES = $(cOBJ) $(cOBJF90) $(cOBJERR)
CutMcK_FFLAGS =  -O3 -march=native -ftree-vectorize -fexpensive-optimizations -ffast-math
CutMcK_FCFLAGS = $(CutMcK_FFLAGS)
CutMcK_LDFLAGS = -lm # -llapack -lblas

# die Includedatei crerr.h muss neu generiert werden, wenn error.txt sich ändert
get_error.f90:	error.txt
		./../make_crerr.sh; 

# zusätzliche clean-rule, um Fortranmodule zu löschen, falls sie existieren
clean-local:
		if [ `ls -1  | sed '/\.mod/!d' | wc -l` -gt 0 ]; then  rm *.mod; fi

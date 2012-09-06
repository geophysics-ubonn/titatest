MODULE modelmod
!!$c 'model.fin'
!!$       
!!$c Andreas Kemna                                            10-Jul-1993
!!$c                                       Letzte Aenderung   24-Oct-1996
!!$
!!$c.....................................................................
!!$c Anzahl der Modellparameter
  INTEGER(KIND = 4),PUBLIC                          :: manz
        
!!$c Zeiger auf Modellparameter
  INTEGER(KIND = 4),PUBLIC,DIMENSION(:),ALLOCATABLE :: mnr

!!!$ >> RM ref model regu
!!$c reference field of reference vector data with some switches
!!!$ which can be enabled in the binary way (tested via BTEST) and thus
!!!$ are additive!
!!!$ wref = 0        -> not use data at all
!!!$ wref = 1        -> use only magnitude data for CRI 
!!!$                    (work only on the REAL part of parameter)
!!!$ wref = wref + 2 -> use magnitude and phase data for CRI and FPI
!!!$ wref = wref + 4 -> reset wref(i) to zero for FPI
  INTEGER(KIND = 4),PUBLIC,DIMENSION(:),ALLOCATABLE :: wref
!!!$ << RM ref model regu
  
END MODULE modelmod

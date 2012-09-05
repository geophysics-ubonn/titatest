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
!!$c reference field of reference vector data (.FALSE. -> exclude, .TRUE. include)
  LOGICAL(KIND = 4),PUBLIC,DIMENSION(:),ALLOCATABLE :: wref
!!!$ << RM ref model regu
  
END MODULE modelmod

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
  
END MODULE modelmod

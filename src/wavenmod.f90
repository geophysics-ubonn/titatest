MODULE wavenmod
!!$
!!$c 'waven.fin'
!!$
!!$c Andreas Kemna                                            07-Oct-1993
!!$c                                       Letzte Aenderung   07-Nov-1997
!!$
!!$c.....................................................................
use alloci
!!!$Anzahl der Wellenzahlwerte
    INTEGER(KIND = 4),PUBLIC                        :: kwnanz

!!!$Schalter steuert Art der Ruecktransformation
!!!$( = 0 : keine Trafo (2D Fall),
!!!$  = 1 : Gauss/Laguerre-Integration )
    INTEGER(KIND = 4),PUBLIC                        :: swrtr

!!!$Wellenzahlwerte
    REAL(prec),PUBLIC,DIMENSION(:),ALLOCATABLE :: kwn

!!!$Wichtungsfaktoren fuer Ruecktransformation
    REAL(prec),PUBLIC,DIMENSION(:),ALLOCATABLE :: kwnwi

!!!$Fuer Ruecktransformation relevanter Abstandsbereich
    REAL(prec),PUBLIC                          :: amin,amax

END MODULE wavenmod

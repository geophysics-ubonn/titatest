MODULE sigmamod
!!$
!!$c 'sigma.fin'
!!$
!!$c Andreas Kemna                                            10-Jul-1993
!!$c                                       Letzte Aenderung   20-Oct-1997
!!$  
!!$c.....................................................................
use alloci
!!$c Referenzleitfaehigkeit
   COMPLEX(prec),PUBLIC                          :: sigma0

!!$c Leitfaehigkeit der Elemente
   COMPLEX(prec),PUBLIC,ALLOCATABLE,DIMENSION(:) :: sigma
!!$c save variable
   COMPLEX(prec),PUBLIC,ALLOCATABLE,DIMENSION(:) :: sigma2

!!$c Background-Werte
   REAL(prec),PUBLIC                             :: bet0,pha0

!!$c Schalter ob "background" Werte eingelesen werden sollen
   LOGICAL,PUBLIC                                     :: lrho0
        
!!$c Schalter ob "starting model" eingelesen werden soll
   LOGICAL,PUBLIC                                     :: lstart
!!$c initial seed fuers verrauschen vom startmodell
   INTEGER(KIND = 4),PUBLIC                           :: iseedpri 
!!$c Rauschen fürs startmodell
   REAL(prec),PUBLIC                             :: modl_stdn


END MODULE sigmamod

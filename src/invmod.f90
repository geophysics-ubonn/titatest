MODULE invmod
!!$c  'inv.fin'
!!$c
!!$c Original version by 
!!$c Andreas Kemna                                            02-Apr-1994
!!$c
!!$c Modified F90 module by Roland Martin
!!$c                                       Letzte Aenderung   27-07-2010
!!$c 
!!$c.....................................................................
use alloci
!!!$ Datenvektor
  COMPLEX(prec),PUBLIC,DIMENSION(:),ALLOCATABLE   :: dat
!!$c Parametervektor
  COMPLEX(prec),PUBLIC,DIMENSION(:),ALLOCATABLE   :: par
!!$c Verbesserungsvektor
  COMPLEX(prec),PUBLIC,DIMENSION(:),ALLOCATABLE   :: dpar
!!$ storage
  COMPLEX(prec),PUBLIC,DIMENSION(:),ALLOCATABLE   :: dpar2
!!$c data vector for difference inversion
  COMPLEX(prec),PUBLIC,DIMENSION(:),ALLOCATABLE   :: d0
!!$c parameter vector for diff inv
  COMPLEX(prec),PUBLIC,DIMENSION(:),ALLOCATABLE   :: m0
!!$c model vector for diff inv
  COMPLEX(prec),PUBLIC,DIMENSION(:),ALLOCATABLE   :: fm0
!!!$ >> RM ref model regu
!!$c parameter vector for reference model regularization
  COMPLEX(prec),PUBLIC,DIMENSION(:),ALLOCATABLE   :: m_ref
!!!$ << RM ref model regu
!!!$ the reference vector is of the same length as par or m0
!!$c data weighting (data covariance) matrix (diag)
  REAL(prec),PUBLIC,DIMENSION(:),ALLOCATABLE      :: wmatd
!!$c storage
  REAL(prec),PUBLIC,DIMENSION(:),ALLOCATABLE      :: wmatd2
!!$c parameter variance
  REAL(prec),PUBLIC                               :: par_vari
!!$c Hilfsfeld, bestimmt welche Daten beruecksichtigt (=1) bzw. nicht
!!$c beruecksichtigt (=0) werden
  INTEGER(KIND = 4),PUBLIC,DIMENSION(:),ALLOCATABLE    :: wdfak
!!$c Schalter ob reine Phaseninversion durchgefuehrt werden soll
  LOGICAL(KIND = 4),PUBLIC                             :: lfpi
!!!$ Auxiliary field accounts for if the background value is used as
!!$ hard reference regularization, i.e. RTR m + I(m - m0), where 
!!!$ each entry of m0 is multiplied with wmfak to cancel values out
  INTEGER (KIND = 4),PUBLIC,DIMENSION(:),ALLOCATABLE   :: wmfak

END MODULE invmod

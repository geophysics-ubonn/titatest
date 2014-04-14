MODULE alloci
! Andreas Kemna                                            24-Jan-1997
!                                       Letzte Aenderung   13-Nov-1997


! Set precision (dp=real64 or qp=real128)
use iso_fortran_env
integer,parameter,public :: prec = Kind(0D0)

! COMPLEX CASE

! Gesamtsteifigkeitsmatrix
  COMPLEX (prec), DIMENSION(:), ALLOCATABLE, PUBLIC     :: a
  COMPLEX (prec), DIMENSION(:,:), ALLOCATABLE, PUBLIC     :: a_mat,b_mat,&
      a_mat_band,a_mat_band_elec
   
  INTEGER,dimension(:),allocatable, public ::          ipiv
!!$Potentialwerte aller Elektrodenlokationen der einzelnen Wellenzahlen
!!$(werden bei der Berechnung der Sensitivitaeten benoetigt)
  COMPLEX (prec), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: kpot
!!$Potentialwerte aller Elektrodenlokationen nach Ruecktransformation
  COMPLEX (prec), DIMENSION(:,:), ALLOCATABLE, PUBLIC   :: hpot
!!$Sensitivitaeten
  COMPLEX (prec), DIMENSION(:,:), ALLOCATABLE, PUBLIC   :: sens
!!$Coverage
  REAL (prec), DIMENSION(:), ALLOCATABLE, PUBLIC   :: csens
!!$DC-CASE
!!$Gesamtsteifigkeitsmatrix
  REAL (prec), DIMENSION(:), ALLOCATABLE, PUBLIC        :: adc
!!$Potentialwerte aller Elektrodenlokationen der einzelnen Wellenzahlen
!!$(werden bei der Berechnung der Sensitivitaeten benoetigt)
  REAL (prec), DIMENSION(:,:,:),ALLOCATABLE, PUBLIC     :: kpotdc
!!$Potentialwerte aller Elektrodenlokationen nach Ruecktransformation
  REAL (prec), DIMENSION(:,:),ALLOCATABLE, PUBLIC       :: hpotdc
!!$Sensitivitaeten
  REAL (prec), DIMENSION(:,:),ALLOCATABLE, PUBLIC       :: sensdc
!!$real symmetric data covariance
  REAL (prec), DIMENSION(:,:), ALLOCATABLE, PUBLIC      :: cov_d
!!$Regularisierungsmatrix
  REAL (prec), DIMENSION(:,:),ALLOCATABLE, PUBLIC       :: smatm
!!$ PSR felder fuer widerstand (r) und phase (p)
  REAL (prec), DIMENSION(:), ALLOCATABLE, PUBLIC        :: rnd_r,rnd_p
!!$real symmetric matrix to compute general inverse
  REAL (prec), DIMENSION(:,:), ALLOCATABLE, PUBLIC   :: ata
!!$general symmetric transpose matrix (regularized) to compute general inverse
  REAL (prec), DIMENSION(:,:), ALLOCATABLE, PUBLIC   :: ata_reg
!!$inverse matrix (may be resolution matrix or the MCM)
  REAL (prec), DIMENSION(:,:), ALLOCATABLE, PUBLIC   :: cov_m
END MODULE alloci

MODULE alloci
!!$Andreas Kemna                                            24-Jan-1997
!!$                                      Letzte Aenderung   13-Nov-1997
!!$COMPLEX CASE
!!$Gesamtsteifigkeitsmatrix
  COMPLEX (KIND(0D0)), DIMENSION(:), ALLOCATABLE, PUBLIC :: a
!!$Potentialwerte aller Elektrodenlokationen der einzelnen Wellenzahlen
!!$(werden bei der Berechnung der Sensitivitaeten benoetigt)
  COMPLEX (KIND(0D0)), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: kpot
!!$Potentialwerte aller Elektrodenlokationen nach Ruecktransformation
  COMPLEX (KIND(0D0)), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: hpot
!!$Sensitivitaeten
  COMPLEX (KIND(0D0)), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: sens
!!$general symmetric transpose matrix to compute general inverse
  COMPLEX (KIND(0D0)), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: ata
!!$general symmetric transpose matrix (regularized) to compute general inverse
  COMPLEX (KIND(0D0)), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: ata_reg
!!$inverse matrix (may be resolution matrix or the MCM)
  COMPLEX (KIND(0D0)), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: cov_m
!!$DC-CASE
!!$Gesamtsteifigkeitsmatrix
  REAL (KIND(0D0)), DIMENSION(:), ALLOCATABLE, PUBLIC :: adc
!!$Potentialwerte aller Elektrodenlokationen der einzelnen Wellenzahlen
!!$(werden bei der Berechnung der Sensitivitaeten benoetigt)
  REAL (KIND(0D0)), DIMENSION(:,:,:),ALLOCATABLE, PUBLIC :: kpotdc
!!$Potentialwerte aller Elektrodenlokationen nach Ruecktransformation
  REAL (KIND(0D0)), DIMENSION(:,:),ALLOCATABLE, PUBLIC :: hpotdc
!!$Sensitivitaeten
  REAL (KIND(0D0)), DIMENSION(:,:),ALLOCATABLE, PUBLIC :: sensdc
!!$real symmetric transpose matrix to compute general inverse
  REAL (KIND(0D0)), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: ata_dc
!!$real symmetric transpose matrix (regularized) to compute general inverse
  REAL (KIND(0D0)), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: ata_reg_dc
!!$real inverse matrix (may be resolution matrix or the MCM)
  REAL (KIND(0D0)), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: cov_m_dc
!!$real symmetric data covariance
  REAL (KIND(0D0)), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: cov_d
!!$Regularisierungsmatrix
  REAL (KIND(0D0)), DIMENSION(:,:),ALLOCATABLE, PUBLIC :: smatm
!!$ PSR felder fuer widerstand (r) und phase (p)
  REAL (KIND(0D0)), DIMENSION(:), ALLOCATABLE, PUBLIC :: rnd_r,rnd_p
END MODULE alloci

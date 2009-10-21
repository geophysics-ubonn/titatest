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
!!$D!!$CASE
!!$Gesamtsteifigkeitsmatrix
  REAL (KIND(0D0)), DIMENSION(:), ALLOCATABLE, PUBLIC :: adc
!!$Potentialwerte aller Elektrodenlokationen der einzelnen Wellenzahlen
!!$(werden bei der Berechnung der Sensitivitaeten benoetigt)
  REAL (KIND(0D0)), DIMENSION(:,:,:),ALLOCATABLE, PUBLIC :: kpotdc
!!$Potentialwerte aller Elektrodenlokationen nach Ruecktransformation
  REAL (KIND(0D0)), DIMENSION(:,:),ALLOCATABLE, PUBLIC :: hpotdc
!!$Sensitivitaeten
  REAL (KIND(0D0)), DIMENSION(:,:),ALLOCATABLE, PUBLIC :: sensdc
!!$Regularisierungsmatrix
  REAL (KIND(0D0)), DIMENSION(:,:),ALLOCATABLE, PUBLIC :: smatm
!!$ PSR felder fuer widerstand (r) und phase (p)
  REAL (KIND(0D0)), DIMENSION(:), ALLOCATABLE, PUBLIC :: rnd_r,rnd_p
END MODULE alloci

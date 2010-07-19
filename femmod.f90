!!$ $Id: fem.f90 1.4 2010/07/19 14:05:36 Roland Martin Exp $
MODULE femmod
!!$ ---------------------------------------------------------
!!$ this module contains some bigger arrays from fem.fin and should 
!!$ replace fem.fin   -.-
!!$ Copyright by Andreas Kemna 2010
!!$ Written by R. Martin 2010
!!$
!!$ included from FEM.FIN!!! COMPLEXs
!!$
!!$ Berechnete Potentialwerte (bzw. Loesungsverktor)
  COMPLEX (KIND(0D0)), DIMENSION(:), ALLOCATABLE, PUBLIC     :: pot
!!$ Analytische berechnete Potentialwerte 
  COMPLEX (KIND(0D0)), DIMENSION(:), ALLOCATABLE, PUBLIC     :: pota
!!$ Konstanten-(bzw. Strom-) Vektor
  COMPLEX (KIND(0D0)), DIMENSION(:), ALLOCATABLE, PUBLIC     :: b
!!$ included from FEM.FIN!!! REALs
  REAL (KIND(0D0)), DIMENSION(:), ALLOCATABLE, PUBLIC        :: bdc
!!$ Skalirerungsfaktor
  REAL (KIND(0D0)), DIMENSION(:), ALLOCATABLE, PUBLIC        :: fak
!!$ Elementbeitraege
  REAL (KIND(0D0)), DIMENSION(:,:,:),ALLOCATABLE, PUBLIC     :: elbg
!!$ Randelementbeitraege
  REAL (KIND(0D0)), DIMENSION(:,:),ALLOCATABLE, PUBLIC       :: relbg
!!$ Konfigurationsfaktoren zur Berechnung der gemischten RB
  REAL (KIND(0D0)), DIMENSION(:,:,:),ALLOCATABLE, PUBLIC     :: kg
!!$ Schalter ob "Gemischte Randbedingung" - Geometrie
  LOGICAL, PUBLIC                                            :: lbeta
!!$ Schalter ob Dirichletsche Randbedingung vorkommt
  LOGICAL, PUBLIC                                            :: lrandb
!!$ Schalter ob "singularity removal" durchgefuehrt werden soll
  LOGICAL, PUBLIC                                            :: lsr
!!$ Schalter ob reine Betragsinversion ("DC") durchgefuehrt werden soll
  LOGICAL, PUBLIC                                            :: ldc
END MODULE femmod

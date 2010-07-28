MODULE datmod
!!$c 'dat.fin'
!!$c
!!$c Original version by 
!!$c Andreas Kemna                                            07-Oct-1993
!!$c
!!$c Modified F90 module by Roland Martin
!!$c                                       Letzte Aenderung   27-07-2010
!!$c.....................................................................
!!$c Anzahl der Messwerte
  INTEGER (KIND = 4),PUBLIC                             :: nanz
!!$c Nummern der Stromelektroden
  INTEGER (KIND = 4), DIMENSION(:),ALLOCATABLE, PUBLIC  :: strnr
!!$c Stromwerte
  REAL (KIND(0D0)), DIMENSION(:),ALLOCATABLE, PUBLIC    :: strom
!!$c Nummern der Spannungselektroden
  INTEGER (KIND = 4), DIMENSION(:),ALLOCATABLE, PUBLIC  :: vnr
!!$c Spannungswerte
  COMPLEX (KIND(0D0)), DIMENSION(:),ALLOCATABLE, PUBLIC :: volt
!!$c Scheinbare Leitfaehigkeiten
  COMPLEX (KIND(0D0)), DIMENSION(:),ALLOCATABLE, PUBLIC :: sigmaa
!!$c storage
  COMPLEX (KIND(0D0)), DIMENSION(:),ALLOCATABLE, PUBLIC :: sgmaa2
!!$c Konfigurationsfaktoren
  REAL (KIND(0D0)), DIMENSION(:),ALLOCATABLE, PUBLIC    :: kfak
!!$c Wichtungsvektor der Phasen
  REAL (KIND(0D0)), DIMENSION(:),ALLOCATABLE, PUBLIC    :: wmatdp
!!$c Wichtungsvektor der Widerstaende
  REAL (KIND(0D0)), DIMENSION(:),ALLOCATABLE, PUBLIC    :: wmatdr
!!$c parameters of resistance error model (dR=stabw0*R+stabm0)
  REAL (KIND(0D0)), PUBLIC                              :: stabw0,stabm0
!!$c parameters of phase error model
!!$c (dp=stabA1*R^stabpB+stabpA2*|P|+stabpA3)
  REAL (KIND(0D0)), PUBLIC                              :: stabp0,stabpA1
  REAL (KIND(0D0)), PUBLIC                              :: stabpB,stabpA2
!!$c Schalter ob 'individual error' (.true.) oder 'uniform weighting'
!!$c (.false.)
  LOGICAL, PUBLIC                                       :: lindiv
!!$c Schalter ob ratio-dataset (.true.)
  LOGICAL, PUBLIC                                       :: lratio
!!$c Schalter ob Polaritaeten gecheckt werden sollen (.true.)
  LOGICAL, PUBLIC                                       :: lpol
!!$c Daten Verrauschen? lnse2 entkoppelt Rauschen und Fehlermodell
  LOGICAL, PUBLIC                                       :: lnse,lnse2
!!$c Initialer Seed der Pseudo Random Sequence (PRS)
  INTEGER (KIND = 4), PUBLIC                            :: iseed
!!$c Noise error model ...
!!$c ... of resistance error model (dnR=nstabw0*nR+nstabm0)
  REAL (KIND(0D0)), PUBLIC                              :: nstabw0,nstabm0
!!$c ... of phase error model 
!!$c (dnp=nstabA1*dnR^nstabpB+nstabpA2*nP+nstabp0)
  REAL (KIND(0D0)), PUBLIC                              :: nstabpB,nstabpA1
  REAL (KIND(0D0)), PUBLIC                              :: nstabpA2,nstabp0
!!$c Anzahl der nicht beruecksichtigten Messwerte
  INTEGER (KIND = 4), PUBLIC                            :: npol

END MODULE datmod

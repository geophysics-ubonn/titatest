        module alloci

c Andreas Kemna                                            24-Jan-1997
c                                       Letzte Aenderung   13-Nov-1997
 
c.....................................................................

c COMPLEX CASE
c Gesamtsteifigkeitsmatrix
        complex         * 16    , dimension(:)
     1                          , ALLOCATABLE, save :: a

c Potentialwerte aller Elektrodenlokationen der einzelnen Wellenzahlen
c (werden bei der Berechnung der Sensitivitaeten benoetigt)
        complex         * 16    , dimension(:,:,:)
     1                          , ALLOCATABLE, save :: kpot

c Potentialwerte aller Elektrodenlokationen nach Ruecktransformation
        complex         * 16    , dimension(:,:)
     1                          , ALLOCATABLE, save :: hpot

c Sensitivitaeten
        complex         * 16    , dimension(:,:)
     1                          , ALLOCATABLE, save :: sens

c DC CASE
c Gesamtsteifigkeitsmatrix
        real            * 8     , dimension(:)
     1                          , ALLOCATABLE, save :: adc

c Potentialwerte aller Elektrodenlokationen der einzelnen Wellenzahlen
c (werden bei der Berechnung der Sensitivitaeten benoetigt)
        real            * 8     , dimension(:,:,:)
     1                          , ALLOCATABLE, save :: kpotdc

c Potentialwerte aller Elektrodenlokationen nach Ruecktransformation
        real            * 8     , dimension(:,:)
     1                          , ALLOCATABLE, save :: hpotdc

c Sensitivitaeten
        real            * 8     , dimension(:,:)
     1                          , ALLOCATABLE, save :: sensdc

        end

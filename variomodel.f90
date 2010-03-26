MODULE variomodel


  IMPLICIT none
  PUBLIC :: rvario
  PUBLIC :: mvario
  PUBLIC :: gvario

  INTEGER(KIND = 4),PRIVATE,SAVE    :: c1 ! set on first call of rvario
  REAL,PRIVATE,SAVE                 :: ome  ! power model exponent
  CHARACTER (30),PRIVATE,SAVE       :: csz  ! string of model type
  REAL,PRIVATE,SAVE                 :: Ix,Iy ! correlation length
CONTAINS

  SUBROUTINE rvario (type,ax,ay,esp_mit,esp_med) 

    INTEGER,INTENT(IN)          :: type ! type of variogram model
    REAL(KIND(0D0)),INTENT(IN)  :: ax,ay ! ax/ay anisotropy coefficients
    REAL(KIND(0D0)),INTENT(IN)  :: esp_mit,esp_med ! aus (besp_elem.for)
    INTEGER                     :: i

    ome = 1.
    c1 = type
    
    IF (ax == 0.) THEN
       Ix=esp_mit
       Iy=esp_mit
       PRINT*,'Choosing mean ESP distance as scale length:',Ix
    ELSE IF (ay == 0.) THEN
       Ix=esp_med
       Iy=esp_med
       PRINT*,'Choosing median ESP distance as scale length:',Ix
    ELSE
       Ix=ax
       Iy=ay
    END IF

    SELECT CASE (c1)
    CASE (1) !spherical
       WRITE (csz,'(a)')'va*(1.5(h/a)-.5*(h/a)**3)'
    CASE (2) ! Gaussian
       WRITE (csz,'(a)')'va*(1-EXP(-(3h/a)**2))'
       Ix = Ix / 3.
       IY = Iy / 3.
    CASE (3) ! power
       PRINT*,'Change power model exponent?[',ome,']'
       READ (*,'(a)')csz
       IF (csz /= '')READ (csz,*)ome
       WRITE (csz,'(a,F2.1)')'va*(h/a)**',ome
    CASE DEFAULT! exponential (default)
       WRITE (csz,'(a)')'va*(1-EXP(-(3h/a)))'
       Ix = Ix / 3.
       IY = Iy / 3.
    END SELECT
    
  END SUBROUTINE rvario

  SUBROUTINE gvario (ax,ay,cvario)
    REAL (KIND(0D0)),INTENT (INOUT) :: ax,ay
    CHARACTER (*)                   :: cvario
    ax = Ix
    ay = Iy
    cvario =  TRIM(csz)
  END SUBROUTINE gvario

  REAL (KIND (0D0)) FUNCTION mvario (lag,varianz)
    ! lag = distance/korrelation (lag)
    REAL (KIND (0D0)),INTENT (IN) :: lag,varianz

    mvario = 0.

    SELECT CASE (c1)
    CASE (1)
       mvario = varianz*(1.5*lag - .5*lag**3.)
       IF (lag > (Ix+Iy)*.5) mvario = varianz
    CASE (2)
       mvario = varianz*(1. - EXP(-lag**2.))
    CASE (3)
       mvario = varianz*lag**ome
    CASE DEFAULT
       mvario = varianz*(1. - EXP(-lag))
!       mvario = varianz*(EXP(-lag))
    END SELECT
    
  END FUNCTION mvario

  REAL (KIND (0D0)) FUNCTION mcova (lag,varianz)
    ! lag = distance/korrelation (lag)
    REAL (KIND (0D0)),INTENT (IN) :: lag,varianz

    SELECT CASE (c1)
    CASE (1)
       mcova = varianz*(1.5*lag - .5*lag**3.)
       IF (lag > (Ix+Iy)*.5) mcova = varianz
    CASE (2)
       mcova = varianz*(1 - EXP(-lag**2.))
    CASE (3)
       mcova = varianz*lag**ome
    CASE DEFAULT
       mcova = varianz*(1 - EXP(-lag))
    END SELECT

    mcova = EXP(-mcova)

  END FUNCTION mcova
  
END MODULE variomodel

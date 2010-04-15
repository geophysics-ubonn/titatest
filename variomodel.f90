MODULE variomodel


  IMPLICIT none
  PUBLIC :: set_vario
  PUBLIC :: get_vario
  PUBLIC :: mvario
  PUBLIC :: mcova

  INTEGER(KIND = 4),PRIVATE,SAVE    :: c1,c2 ! is set on first call
!!!$ c1 accounts for the variogram model
!!!$ c2 accounts for the covariance model
  REAL,PRIVATE,SAVE                 :: omev,omec  
!!$! power model exponent for variogram and covariance
  CHARACTER (30),PRIVATE,SAVE       :: cszv,cszc  
!!$ strings of model type, cszv for variogram model
!!$ cszc for covariance model (can be decoupeled)
  REAL,PRIVATE,SAVE                 :: Ix,Iy
!!!$ this may look a littel strange but in fact there are some
!!!$ variogram models (Exponential and Gauss) which need 
!!!$ the Correlation length to be /3 of it. so be careful to mistake
!!!$ these numbers wrong.--

  PRIVATE 
CONTAINS

  SUBROUTINE set_vario (type,ax,ay,esp_mit,esp_med) 

    INTEGER,INTENT(IN)          :: type ! type of variogram model
    REAL(KIND(0D0)),INTENT(IN)  :: ax,ay ! ax/ay anisotropy coefficients
    REAL(KIND(0D0)),INTENT(IN)  :: esp_mit,esp_med ! aus (besp_elem.for)
    INTEGER                     :: i

    omev = 1.5 ! 0<ome<2
    omec = omev
    
    c2 = INT(type/10)
    c1 = type-c2*10

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
       WRITE (cszv,'(a)')'va*(1-EXP(-(3h/a)**2))'
    CASE (2) ! Gaussian
       PRINT*,'Change power model exponent?[',omev,']'
       READ (*,'(a)')cszv
       IF (cszv /= '')READ (cszv,*)omev
       WRITE (cszv,'(a,F2.1)')'va*(h/a)**',omev
    CASE (3) ! power
       WRITE (cszv,'(a)')'va*(1.5(h/a)-.5*(h/a)**3),va'
    CASE DEFAULT! exponential (default)
       WRITE (cszv,'(a)')'va*(1-EXP(-(3h/a))) (default)'
    END SELECT

    SELECT CASE (c2)
    CASE (1) !Continous Exponential
       WRITE (cszc,'(a)')'va*EXP(-3h/a)'
    CASE (2) !Continous Gaussian
       WRITE (cszc,'(a)')'va*EXP(-(3h/a)**2)'
    CASE (3) ! power
       PRINT*,'Change power model exponent?[',omec,']'
       READ (*,'(a)')cszc
       IF (cszc /= '')READ (cszc,*)omec
       WRITE (cszc,'(a,F2.1,a)')'EXP(-va*(h/a)**',omec,')'
    CASE (4)
       WRITE (cszc,'(a)')'EXP(-va*(1.5(h/a)-.5*(h/a)**3),-va)'
    CASE DEFAULT! exponential (default)
       WRITE (cszc,'(a)')'EXP(-variogram(h)) (default)'
    END SELECT
    
  END SUBROUTINE set_vario

  SUBROUTINE get_vario (ax,ay,csz,type)
    INTEGER,INTENT(IN)              :: type 
!!$! which info type=0 -> variogram type = 1->covariance
    REAL (KIND(0D0)),INTENT (INOUT) :: ax,ay
    CHARACTER (*)                   :: csz

    ax = Ix
    ay = Iy

    SELECT CASE (type)
    CASE (0)
       csz =  TRIM(cszv)
    CASE (1)
       csz =  TRIM(cszc)
    END SELECT

  END SUBROUTINE get_vario

  REAL (KIND (0D0)) FUNCTION mvario (lagx,lagy,varianz)
    ! lag = distance/korrelation (lag)
    REAL (KIND (0D0)),INTENT (IN) :: lagx,lagy,varianz
    REAL (KIND (0D0))             :: dum,r

    mvario = 0.
    dum = SQRT((lagx / Ix)**2. + (lagy / Iy)**2.)

    SELECT CASE (c1)
    CASE (1)
       dum = SQRT((lagx * 3. / Ix)**2. + (lagy * 3. / Iy)**2.)
       mvario = varianz*(1. - EXP(-dum**2.))
    CASE (2)
       mvario = varianz*dum**omev
    CASE (3)
       r = SQRT(lagx**2. + lagy**2.)
       IF (r < (Ix + Iy)*.5) THEN
          mvario = varianz*(1.5*dum - .5*dum**3.)
       ELSE
          mvario = varianz
       END IF
    CASE DEFAULT
       dum = SQRT((lagx * 3. / Ix)**2. + (lagy * 3. / Iy)**2.)
       mvario = varianz*(1. - EXP(-dum))
    END SELECT
    
  END FUNCTION mvario

  REAL (KIND (0D0)) FUNCTION mcova (lagx,lagy,varianz)
    ! lag = distance/korrelation (lag)
    REAL (KIND (0D0)),INTENT (IN) :: lagx,lagy,varianz
    REAL (KIND (0D0))             :: dum,r
    
    dum = SQRT((lagx / Ix)**2. + (lagy / Iy)**2.)

    SELECT CASE (c2)

    CASE (1)
       dum = SQRT((lagx * 3. / Ix)**2. + (lagy * 3. / Iy)**2.)
       mcova = varianz*EXP(-dum)
    CASE (2)
       dum = SQRT((lagx * 3. / Ix)**2. + (lagy * 3. / Iy)**2.)
       mcova = varianz*EXP(-dum**2.)
    CASE (3)
       mcova = varianz*dum**omec
       mcova = EXP(-mcova)
    CASE (4)
       r = SQRT(lagx**2. + lagy**2.)
       IF (r < (Ix + Iy)*.5) THEN
          mcova = varianz*(1.5*dum - .5*dum**3.)
       ELSE
          mcova = varianz
       END IF
       mcova = EXP(-mcova)
    CASE DEFAULT
       mcova = EXP(-mvario(lagx,lagy,varianz))
    END SELECT

  END FUNCTION mcova
  
END MODULE variomodel

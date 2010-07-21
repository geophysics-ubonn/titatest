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
  REAL,PRIVATE,SAVE                 :: axs,ays
!!!$ this may look a littel strange but in fact there are some
!!!$ variogram models (Exponential and Gauss) which need 
!!!$ the Correlation length to be /3 of it. so be careful to mistake
!!!$ these numbers wrong.--

  PRIVATE 
CONTAINS

  SUBROUTINE set_vario (type,ax,ay,esp_mit,esp_med) 
! is called at the beginning
    INTEGER,INTENT(IN)          :: type ! type of variogram model
    REAL(KIND(0D0)),INTENT(IN)  :: ax,ay ! ax/ay anisotropy coefficients
    REAL(KIND(0D0)),INTENT(IN)  :: esp_mit,esp_med ! aus (besp_elem.for)
    INTEGER                     :: i

    omev = 1.5 ! 0<ome<2
    omec = omev
    
    c2 = INT(type/10)
    c1 = type-c2*10

    IF (ax == 0.) THEN ! taking default values if no value
       Ix = esp_mit ! arithmetical mean
       Iy = esp_mit
       PRINT*,'Choosing mean ESP distance as scale length:',Ix
    ELSE IF (ay == 0.) THEN
       Ix = esp_med ! median
       Iy = esp_med
       PRINT*,'Choosing median ESP distance as scale length:',Ix
    ELSE
       Ix = ax
       axs = ax
       Iy = ay
       ays = ay
    END IF

    SELECT CASE (c1) ! string for variogram function
    CASE (1) !Gaussian 
       Ix = Ix/9.
       Iy = Iy/9.
       WRITE (cszv,'(a)')'va*(1-EXP(-(3h/a)**2))'
    CASE (2) ! Spherical
       WRITE (cszv,'(a)')'va((1.5(h/a)-.5*(h/a)**3),1)'
    CASE (3) ! power
       READ (*,'(a)')cszv
       IF (cszv /= '')READ (cszv,*)omev
       WRITE (cszv,'(a,F2.1)')'va*(h/a)**',omev
    CASE DEFAULT! exponential (default)
       Ix = Ix/3.
       Iy = Iy/3.
       WRITE (cszv,'(a)')'va*(1-EXP(-(3h/a))) (default)'
    END SELECT

    SELECT CASE (c2)
    CASE (1) !Gaussian
       WRITE (cszc,'(a)')'va*EXP(-(3h/a)**2)'
    CASE (2) !Spherical
       WRITE (cszc,'(a)')'(va*(1-1.5(h/a)+.5(h/a)**3),0)'
    CASE (3) !Power
       PRINT*,'Change power model exponent?[',omec,']'
       READ (*,'(a)')cszc
       IF (cszc /= '')READ (cszc,*)omec
       WRITE (cszc,'(a,F2.1,a)')'EXP(-va*(h/a)**',omec,')'
    CASE (4)!Lemma
       WRITE (cszc,'(a)')'EXP(-variogram(h))'
    CASE DEFAULT!Exponential1 (default)
       WRITE (cszc,'(a)')'va*EXP(-3h/a) (default)'
    END SELECT
    
  END SUBROUTINE set_vario

  SUBROUTINE get_vario (ax,ay,csz,type)
    INTEGER,INTENT(IN)              :: type 
!!$! which info type=0 -> variogram type = 1->covariance
    REAL (KIND(0D0)),INTENT (INOUT) :: ax,ay
    CHARACTER (*)                   :: csz
! gives back the integral scale used for the variogram
    ax = axs
    ay = ays

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
    REAL (KIND (0D0))             :: rh,ra,r

    mvario = 0.
    r = SQRT((lagx / Ix)**2. + (lagy / Iy)**2.)
    rh = SQRT(lagx**2. + lagy**2.)
    ra = SQRT(Ix**2. + Iy**2.)

    SELECT CASE (c1)
    CASE (1)
       mvario = varianz * (1. - EXP(-r**2.))
    CASE (2)
       IF (rh <= ra) THEN
          mvario = varianz * (1.5*r - .5*r**3.)
       ELSE
          mvario = varianz
       END IF
    CASE (3)
       mvario = varianz * r**omec
    CASE DEFAULT
       mvario = varianz*(1. - EXP(-r))
    END SELECT
    
  END FUNCTION mvario

  REAL (KIND (0D0)) FUNCTION mcova (lagx,lagy,varianz)
    ! lag = distance/korrelation (lag)
    REAL (KIND (0D0)),INTENT (IN) :: lagx,lagy,varianz
    REAL (KIND (0D0))             :: r,rh,ra ! distances
    
    r = SQRT((lagx / Ix)**2. + (lagy / Iy)**2.) ! not sure of this
    rh = SQRT(lagx**2. + lagy**2.)
    ra = SQRT(Ix**2. + Iy**2.)
    mcova = 0.

    SELECT CASE (c2)

    CASE (1)
       mcova = varianz*EXP(-r**2.)
    CASE (2)
       IF (rh <= ra) THEN
          mcova = varianz * (1. - 1.5*r + .5*r**3.)
       ELSE
          mcova = 0.
       END IF
    CASE (3)
       mcova = varianz*r**omec
       mcova = EXP(-mcova)
    CASE (4)
       mcova = EXP(-mvario(lagx,lagy,varianz))
    CASE DEFAULT
       mcova = varianz * EXP(-r)
    END SELECT

  END FUNCTION mcova
  
END MODULE variomodel

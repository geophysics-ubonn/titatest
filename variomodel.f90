MODULE variomodel


  IMPLICIT none
  PUBLIC :: set_vario
  PUBLIC :: get_vario
  PUBLIC :: mvario
  PUBLIC :: mcova

  INTEGER(KIND = 4),PRIVATE,SAVE    :: c1,c2 ! is set on first call
!!!$ c1 accounts for the variogram model
!!!$ c2 accounts for the covariance model
  REAL(KIND(0D0)),PRIVATE,SAVE      :: omev,omec
!!$! power model exponent for variogram and covariance
  REAL(KIND(0D0)),PRIVATE,SAVE      :: tfac
!!$! exponent factor for covariance function according to exp(-tfac(variogram))
  CHARACTER (30),PRIVATE,SAVE       :: cszv,cszc  
!!$ strings of model type, cszv for variogram model
!!$ cszc for covariance model (can be decoupeled)
  REAL(KIND(0D0)),PRIVATE,SAVE      :: Ix_v,Iy_v
!!$ Integral scales for variogram function
  REAL(KIND(0D0)),PRIVATE,SAVE      :: Ix_c,Iy_c
!!$ Integral scales for covariance function
!!!$ this may look a littel strange but in fact there are some
!!!$ variogram models (Exponential and Gauss) which need 
!!!$ the Correlation length to be /3 of it. so be careful to mistake
!!!$ these numbers wrong.--
  REAL(KIND(0D0)),PRIVATE,SAVE      :: axs,ays
!!$ True integral scale from user..

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
    tfac = omev

    c2 = INT(type/10)
    c1 = type-c2*10

    IF (ax == 0.) THEN ! taking default values if no value
       axs = esp_mit ! arithmetical mean of bnachbar
       ays = esp_mit
       PRINT*,'Choosing mean ESP distance as scale length:',esp_mit
    ELSE IF (ay == 0.) THEN
       axs = esp_med ! median
       ays = esp_med
       PRINT*,'Choosing median ESP distance as scale length:',esp_med
    ELSE
       axs = ax ! save user values
       ays = ay
    END IF

    Ix_v = axs;Ix_c = axs ! sets integral scales--
    Iy_v = ays;Iy_c = ays

    SELECT CASE (c1) ! string for variogram function
    CASE (1) !Gaussian 
       Ix_v = Ix_v/9. ! scale is changed to match GSlib standard
       Iy_v = Iy_v/9.
       WRITE (cszv,'(a)')'va(1-EXP(-(3h/a)**2))'
    CASE (2) ! Spherical
       WRITE (cszv,'(a)')'va((1.5(h/a)-.5(h/a)**3),1)'
    CASE (3) ! power
       READ (*,'(a)')cszv
       IF (cszv /= '')READ (cszv,*)omev
       WRITE (cszv,'(a,F3.1)')'va(h/a)**',omev
    CASE DEFAULT! exponential
       Ix_v = Ix_v/3.! scale is changed to match GSlib standard
       Iy_v = Iy_v/3.
       WRITE (cszv,'(a)')'va(1-EXP(-(3h/a)))'
    END SELECT

    SELECT CASE (c2)
    CASE (1) !Gaussian
       Ix_c = Ix_c/9.! scale is changed to match GSlib standard
       Iy_c = Iy_c/9.
       WRITE (cszc,'(a)')'vaEXP(-(3h/a)**2)'
    CASE (2) !Spherical
       WRITE (cszc,'(a)')'va((1-1.5(h/a)+.5(h/a)**3),0)'
    CASE (3) !Power
       PRINT*,'Change power model exponent?[',omec,']'
       READ (*,'(a)')cszc
       IF (cszc /= '')READ (cszc,*)omec
       WRITE (cszc,'(a,F3.1,a)')'EXP(-va*(h/a)**',omec,')'
    CASE (4)!Lemma
       PRINT*,'Change exponent factor?[',tfac,']'
       READ (*,'(a)')cszc
       IF (cszc /= '')READ (cszc,*)tfac
       WRITE (cszc,'(a,F3.1,a)')'EXP(-',tfac,'*variogram(h))'
    CASE DEFAULT!Exponential1
       Ix_c = Ix_c/3.
       Iy_c = Iy_c/3.
       WRITE (cszc,'(a)')'va*EXP(-3h/a)'
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
    REAL (KIND (0D0))             :: r,r2 ! distances

    mvario = 0.

    r = SQRT((lagx / Ix_v)**2. + (lagy / Iy_v)**2.)
    r2 = r*r

    SELECT CASE (c1)
    CASE (1)
       mvario = varianz * (1. - EXP(-r2))
    CASE (2)
       IF (r < 1.) THEN
          mvario = varianz * (r * (1.5 - .5*r2))
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
    REAL (KIND (0D0))             :: r,r2 ! distances
    
    mcova = 0.
    
    r = SQRT((lagx / Ix_c)**2. + (lagy / Iy_c)**2.) ! not sure of this
    r2 = r**1.999999 ! 9 up to seventh digit 
    !! and C of gauss model will not be pos def
    r2 = r*r ! just to be sure -.-
    
    SELECT CASE (c2)
    CASE (1)
       mcova = varianz * EXP(-r2)
    CASE (2)
       IF ( r < 1. ) THEN
          mcova = varianz * (1. - r * (1.5 - .5*r2) )
       ELSE
          mcova = 0.
       END IF
    CASE (3)
       mcova = varianz*r**omec
       mcova = EXP(-mcova)
    CASE (4)
       mcova = EXP(-tfac*mvario(lagx,lagy,varianz))
    CASE DEFAULT
       mcova = varianz * EXP(-r)
    END SELECT
    !    IF (mcova < 1.d-9) mcova = 0.
  END FUNCTION mcova
  
END MODULE variomodel

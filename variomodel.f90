MODULE variomodel


  IMPLICIT none
  PUBLIC :: set_vario
  PUBLIC :: get_vario
  PUBLIC :: mvario
  PUBLIC :: mcova

!!!$ c1 accounts for the variogram model
!!!$ c2 accounts for the covariance model
  INTEGER(KIND = 4),PRIVATE,SAVE    :: c1,c2 ! switches are set on first call
!!!$ the switches are needed internal to discriminate between the various 
!!!$ variogram and covariance functions.
  REAL(KIND(0D0)),PRIVATE,SAVE      :: omev,omec
!!$! power model exponent for variogram and covariance
  REAL(KIND(0D0)),PRIVATE,SAVE      :: tfac
!!$! exponent factor for covariance function according to exp(-tfac(variogram))
  CHARACTER (30),PRIVATE,SAVE       :: cszv,cszc  
!!$ strings of model type, cszv for variogram model
!!$ cszc for covariance model (can be decoupeled)
  REAL(KIND(0D0)),PRIVATE,SAVE      :: Ix_v,Iy_v
!!$ Integral scales (range) for variogram function
  REAL(KIND(0D0)),PRIVATE,SAVE      :: Ix_c,Iy_c
!!$ Integral scales for covariance function.
!!!$ It may look a littel strange but in fact there are some
!!!$ variogram models (Exponential and Gauss) which need 
!!!$ the Correlation length to be /3 (EXP) or /3^2 (GAU) of it. 
!!!$ So, be careful to mistake these numbers wrong. 
!!!$ See also the GSlib manual
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
    tfac = omev ! just a "educated guess"

    c2 = INT(type/10) ! this sets the switches from outside
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
       axs = ax ! save user values into SAVED
       ays = ay
    END IF

    Ix_v = axs;Ix_c = axs ! sets integral scales (range) for 
    Iy_v = ays;Iy_c = ays ! internal usage

    SELECT CASE (c1) ! string for variogram function
    CASE (1) !Gaussian 
       Ix_v = Ix_v / 3D0 ! scale is changed to match GSlib standard
       Iy_v = Iy_v / 3D0 ! 4/7 is from Kitanidis..
       WRITE (cszv,'(a)')'va(1-EXP(-(3h/a)^2))'
    CASE (2) ! Spherical
       WRITE (cszv,'(a)')'va((1.5(h/a)-.5(h/a)^3),1)'
    CASE (3) ! Power
       PRINT*,'Change power model exponent?[',omev,']'
       READ (*,'(a)')cszv
       IF (cszv /= '')READ (cszv,*)omev
       WRITE (cszv,'(a,F3.1)')'va(h/a)^',omev
    CASE DEFAULT! exponential
       Ix_v = Ix_v / 3D0! scale is changed to match GSlib standard
       Iy_v = Iy_v / 3D0
       WRITE (cszv,'(a)')'va(1-EXP(-(3h/a)))'
    END SELECT

    SELECT CASE (c2)
    CASE (1) !Gaussian
       Ix_c = Ix_c / 3D0! scale is changed to match GSlib standard
       Iy_c = Iy_c / 3D0
       WRITE (cszc,'(a)')'vaEXP(-(3h/a)^2)'
    CASE (2) !Spherical
       WRITE (cszc,'(a)')'va((1-1.5(h/a)+.5(h/a)^3),0)'
    CASE (3) !Power
       PRINT*,'Change power model exponent?[',omec,']'
       READ (*,'(a)')cszc
       IF (cszc /= '')READ (cszc,*)omec
       WRITE (cszc,'(a,F3.1,a)')'EXP(-va*(h/a)^',omec,')'
    CASE (4)!Lemma
       PRINT*,'Change exponent factor?[',tfac,']'
       READ (*,'(a)')cszc
       IF (cszc /= '')READ (cszc,*)tfac
       WRITE (cszc,'(a,F3.1,a)')'EXP(-',tfac,'*variogram(h))'
    CASE DEFAULT!Exponential1
       Ix_c = Ix_c / 3D0
       Iy_c = Iy_c / 3D0
       WRITE (cszc,'(a)')'va*EXP(-3h/a)'
    END SELECT
    
  END SUBROUTINE set_vario

  SUBROUTINE get_vario (ax,ay,csz,type)
    INTEGER,INTENT(IN)              :: type 
!!$! which info type=0 -> variogram type = 1->covariance
    REAL (KIND(0D0)),INTENT (OUT) :: ax,ay
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
    r2 = r*r ! just to be sure to have less numerical issues

    SELECT CASE (c1)
    CASE (1)
       mvario = varianz * (1. - EXP(-r2)) !from GSlib
    CASE (2)
       IF (r < 1.) THEN
          mvario = varianz * (r * (1.5 - .5*r2)) !from GSlib
       ELSE
          mvario = varianz
       END IF
    CASE (3)
       mvario = varianz * r**omec !from GSlib
    CASE DEFAULT
       mvario = varianz*(1. - EXP(-r)) !from GSlib
    END SELECT
    
  END FUNCTION mvario

  REAL (KIND (0D0)) FUNCTION mcova (lagx,lagy,varianz)
    ! lag = distance/korrelation (lag) varianz = variance (sill)
    REAL (KIND (0D0)),INTENT (IN) :: lagx,lagy,varianz
    REAL (KIND (0D0))             :: r,r2 ! distances
    mcova = 0.
    
    r = SQRT((lagx / Ix_c)**2. + (lagy / Iy_c)**2.)
    
    r2 = r*r 

    r2 = r**1.99999 ! this is odd.. if we put it to 2., C_m is no longer
!!$ pos definite.. but with 1.99999, which is in fact nearly 2. it is ok
!!$ happens not all the time but sometimes..
!!$ the digit can vary up to 7 digits behind the dot. With 8 digits it becomes
!!$ somehow unstable..
!!$ I took 5 digits to be sure.
!!$ To me this sounds like a DOUBLE_PRECISION / SINGLE_PRECISION Problem..

    SELECT CASE (c2)
    CASE (1)
       mcova = varianz * EXP(-r2) !from GSlib
    CASE (2)
       IF ( r < 1. ) THEN
          mcova = varianz * (1. - r * (1.5 - .5*r2) ) !from GSlib
       ELSE
          mcova = 0.
       END IF
    CASE (3)
       mcova = varianz*r**omec !from GSlib
       mcova = EXP(-mcova) ! own interpretation
    CASE (4)
       mcova = EXP(-tfac*mvario(lagx,lagy,varianz)) ! this is from a lemma
    CASE DEFAULT
       mcova = varianz * EXP(-r) !from GSlib
    END SELECT

  END FUNCTION mcova
  
END MODULE variomodel

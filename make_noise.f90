!!$ $Id: make_noise.f90 1.4 2008/07/29 14:05:36 Roland Martin Exp $
MODULE Make_noise
!!$------------------------------------------------------------------------
!!$  make_noise.f90
!!$
!!$  Contains subroutines and functions to generate ensembles of
!!$  random numbers which shoud be portable (e.g. Press, et al. 1998)
!!$  the distributions are of popular algorithms
!!$
!!$------------------------------------------------------------------------
  IMPLICIT none
  PUBLIC :: Random_Init
  PUBLIC :: Random_Gauss
  PUBLIC :: Random_Exponential
  PUBLIC :: Random_Breitwigner
  PRIVATE :: Random_Draw
  PRIVATE :: Random_Double
  PRIVATE :: Random_Int
  
!!$ now the most important variables..
  INTEGER(KIND = 4),PARAMETER,PRIVATE:: N_RND=97 
  ! number of pseudo random numbers
  REAL(KIND(0D0)),SAVE,DIMENSION(N_RND),PRIVATE :: Rnd 
  !  pseudo random ensemble
  INTEGER(KIND = 4),SAVE,PRIVATE :: ix1,ix2,ix3 
  ! huge numbers which are stored and altered thoroughly
  INTEGER(KIND = 4),PARAMETER,PRIVATE:: M1=259200
  INTEGER(KIND = 4),PARAMETER,PRIVATE:: M2=134456
  INTEGER(KIND = 4),PARAMETER,PRIVATE:: M3=243000
  INTEGER(KIND = 4),PARAMETER,PRIVATE:: IA1=7141
  INTEGER(KIND = 4),PARAMETER,PRIVATE:: IA2=8121
  INTEGER(KIND = 4),PARAMETER,PRIVATE:: IA3=4561
  INTEGER(KIND = 4),PARAMETER,PRIVATE:: IC1=54773
  INTEGER(KIND = 4),PARAMETER,PRIVATE:: IC2=28441
  INTEGER(KIND = 4),PARAMETER,PRIVATE:: IC3=51349
  ! a bunch of magic numbers...
  REAL(KIND(0D0)),PARAMETER,PRIVATE:: RM1=1./M1
  REAL(KIND(0D0)),PARAMETER,PRIVATE:: RM2=1./M2

  CONTAINS
!!$------------------------------------------------------------------------
!!$ initialize portable pseudo random numbers 
!!$ (set the pseudo random numbers)
!!$
    SUBROUTINE Random_Init(iseed)
      INTEGER ( KIND = 4 ),OPTIONAL :: iseed
      INTEGER ( KIND = 4 ) :: i

      IF (.NOT.PRESENT(iseed)) iseed = 1 !default seed

      ix1=MOD(IC1+iseed,M1)
      ix1=MOD(IA1*ix1+IC1,M1)
      ix2=MOD(ix1,M2)
      ix1=MOD(IA1*ix1+ic1,M1)
      ix3=MOD(ix1,M3)
      DO i=1,N_RND
         ix1=MOD(IA1*ix1+IC1,M1)
         ix2=MOD(IA2*ix2+IC2,M2)
         Rnd(i)=(DBLE(ix1)+DBLE(ix2)*RM2)*RM1
      END DO
    END SUBROUTINE Random_Init
!!$------------------------------------------------------------------------
!!$ Draw a random number from the pseudo random sequence
!!$
    REAL (KIND(0D0)) FUNCTION Random_Draw()
      INTEGER ( KIND = 4 ) :: i

      DO
         ix1=MOD(IA1*ix1+IC1,M1)
         ix2=MOD(IA2*ix2+IC2,M2)
         ix3=MOD(IA3*ix3+IC3,M3)
         i=1+(N_RND*ix3)/M3
         IF (i>N_RND.OR.i<1) CYCLE
         Random_Draw=Rnd(i)
         Rnd(i)=(DBLE(ix1)+DBLE(ix2)*RM2)*RM1
         EXIT
      END DO
    END FUNCTION Random_Draw
!!$------------------------------------------------------------------------
!!$ Flat (uniform) distribution
!!$ Returns a uniformly distributed random real number between [min,max]
!!$
    REAL (KIND (0D0)) FUNCTION Random_Double(min, max)
      REAL (KIND (0D0)), INTENT(IN), OPTIONAL :: min, max
      REAL (KIND (0D0)) :: x
      
      x = Random_Draw()
      
      IF(PRESENT(min) .AND. PRESENT(max)) THEN
         Random_Double = min + ABS(max-min)*x
      ELSE
         Random_Double = x
     END IF
     
   END FUNCTION Random_Double
!!$------------------------------------------------------------------------
!!$ Integer (uniform) distribution
!!$ Returns a uniformly distributed random integer number between [min,max]
!!$
   INTEGER FUNCTION Random_Int(min, max)
   INTEGER, INTENT(IN), OPTIONAL :: min, max
   REAL (KIND (0D0)) :: x

     x = Random_Draw()

     IF(PRESENT(min) .AND. PRESENT(max)) THEN
        Random_Int = min + INT(max*x)
     ELSE
        Random_Int = INT(100*x)        
     END IF
     
   END FUNCTION Random_Int
!!$------------------------------------------------------------------------
!!$ Gaussian Distribution
!!$ Returns a normally distributed deviate with mean and sigma
!!$ The routine uses the Box-Muller transformation of uniform deviates.
!!$
   REAL (KIND (0D0)) FUNCTION Random_Gauss(mean, sigma)
   INTEGER, INTENT(IN), OPTIONAL :: mean, sigma
   REAL (KIND (0D0)) :: x, y, z

     DO
        x = 2.0 * Random_Double() - 1.0
        y = 2.0 * Random_Double() - 1.0
        z = x*x + y*y
        if( z <= 1.0 ) exit
     END DO

     IF(PRESENT(mean) .AND. PRESENT(sigma)) THEN
        Random_Gauss = mean + sigma*x*sqrt(-2.0*log(z)/z)
     ELSE
        Random_Gauss = x*sqrt(-2.0*log(z)/z)
     END IF

   END FUNCTION Random_Gauss
!!$------------------------------------------------------------------------
!!$ Exponential (decay) distribution
!!$ Returns a random number between times t1 and t2 
!!$ according to f(t) = exp (-t/tau)
!!$
   REAL (KIND (0D0)) FUNCTION Random_Exponential(tau, tmin, tmax)
   REAL (KIND (0D0)), INTENT(IN) :: tau
   REAL (KIND (0D0)), INTENT(IN), OPTIONAL :: tmin, tmax
   REAL (KIND (0D0)) :: r1, r2

     IF(PRESENT(tmin) .AND. PRESENT(tmax)) THEN
        r1 =  exp(-tmin/tau)
        r2 =  exp(-tmax/tau)
     ELSE
        r1 = 1.0
        r2 = 0.0
     END IF

     Random_Exponential = -tau*log(r2 + Random_Double() * (r1-r2) )
   
   END FUNCTION Random_Exponential
!!$------------------------------------------------------------------------
!!$ Breit-Wigner Distribution
!!$ Returns a random number from a Breit-Wigner distribution 
!!$ for center mean Full Width Half Maximum fwhm
!!$
   REAL (KIND (0D0)) FUNCTION Random_BreitWigner(mean, fwhm)
   REAL (KIND (0D0)), INTENT(IN), OPTIONAL :: mean, fwhm
   REAL (KIND (0D0)) :: x, y, z
   
     DO
        x = 2.0 * Random_Double() - 1.0
        y = 2.0 * Random_Double() - 1.0
        z = x*x + y*y
        if( z <= 1.0 ) exit
     END DO

     IF(PRESENT(mean) .AND. PRESENT(fwhm)) THEN
        Random_BreitWigner = mean + 0.5*fwhm*x/y
     ELSE
        Random_BreitWigner = 0.5*x/y
     END IF

   END FUNCTION Random_BreitWigner

 END MODULE Make_noise

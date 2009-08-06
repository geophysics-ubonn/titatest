!!$ $Id: noise.f90,v 1.4 2008/07/29 14:05:36 Roland Martin Exp $
MODULE make_noise

  IMPLICIT none

  PUBLIC :: dp_noise
  PUBLIC :: initrand
  PUBLIC :: ran1
  PUBLIC :: ran0

CONTAINS

  REAL (KIND(0D0)) FUNCTION dp_noise(idum)
    LOGICAL,SAVE :: done
    INTEGER ( KIND = 4 ), INTENT (INOUT) :: idum
    REAL (KIND(0D0)) :: v1,v2,r,fact
    REAL (KIND(0D0)), SAVE :: gset

    IF ( .not. done ) THEN
       DO WHILE ( .not. done )
          v1 = 2.*RAN1(idum) - 1.
          v2 = 2.*RAN1(idum) - 1.
          r = v1**2 + v2**2
          IF( r >= 1 ) CYCLE
          fact = SQRT( -2.*LOG( r ) / r )
          gset = v1*fact
          dp_noise = v2*fact
          done = .true.
       END DO
    ELSE
       dp_noise = gset
       done = .FALSE.
    END IF
  END FUNCTION DP_NOISE

  REAL (KIND(0D0)) FUNCTION RAN0(idum)
    INTEGER ( KIND = 4 ) :: idum 
    EXTERNAL :: SRAND
    REAL (KIND(0D0)) :: V(97),Y,DUM,rand
    INTEGER ( KIND = 4 ) :: J,IFF,ISEED
    LOGICAL,SAVE :: done
    done=.false.
    DO WHILE (.NOT.done)
       IF (IDUM.LT.0.OR..NOT.done) THEN
          ISEED=ABS(IDUM)
          CALL RANDOM_SEED(ISEED)
          IDUM=1
          DO J=1,97
             DUM=RAND()
          END DO
          DO J=1,97
             V(J)=RAND()
          END DO
          Y=RAND()
       END IF
       J=1+INT(97.*Y)
       IF (J.GT.97.OR.J.LT.1) CYCLE
       Y=V(J)
       RAN0=Y
       V(J)=RAND()
       done=.not.done
    END DO
  END FUNCTION RAN0

!!$ portable random numbers...
  REAL (KIND(0D0)) FUNCTION RAN1(idum)
    REAL (KIND(0D0)), SAVE, DIMENSION(97) :: r
    INTEGER ( KIND = 4 ),PARAMETER:: M1=259200
    INTEGER ( KIND = 4 ),PARAMETER:: M2=134456
    INTEGER ( KIND = 4 ),PARAMETER:: M3=243000
    INTEGER ( KIND = 4 ),PARAMETER:: IA1=7141
    INTEGER ( KIND = 4 ),PARAMETER:: IA2=8121
    INTEGER ( KIND = 4 ),PARAMETER:: IA3=4561
    INTEGER ( KIND = 4 ),PARAMETER:: IC1=54773
    INTEGER ( KIND = 4 ),PARAMETER:: IC2=28441
    INTEGER ( KIND = 4 ),PARAMETER:: IC3=51349
    REAL (KIND(0D0)),PARAMETER:: RM1=1./M1
    REAL (KIND(0D0)),PARAMETER:: RM2=1./M2
    INTEGER ( KIND = 4 ),INTENT(INOUT) :: idum
    INTEGER ( KIND = 4 )::j
    INTEGER ( KIND = 4 ),SAVE :: ix1,ix2,ix3

    DO 
       IF (idum.LT.0) THEN
          ix1=MOD(IC1-idum,M1)
          ix1=MOD(IA1*ix1+IC1,M1)
          ix2=MOD(ix1,M2)
          ix1=MOD(IA1*ix1+ic1,M1)
          ix3=MOD(ix1,M3)
          DO j=1,97
             ix1=MOD(IA1*ix1+IC1,M1)
             ix2=MOD(IA2*ix2+IC2,M2)
             r(j)=(FLOAT(ix1)+FLOAT(ix2)*RM2)*RM1
          END DO
          idum=1!ist initialisiert...
       END IF
       ix1=MOD(IA1*ix1+IC1,M1)
       ix2=MOD(IA2*ix2+IC2,M2)
       ix3=MOD(IA3*ix3+IC3,M3)
       j=1+(97*ix3)/M3
       IF (j>97.OR.j<1) CYCLE
       ran1=r(j)
       r(j)=(FLOAT(ix1)+FLOAT(ix2)*RM2)*RM1
       EXIT
    END DO
  END FUNCTION RAN1

  INTEGER ( KIND = 4 ) FUNCTION INITRAND()
    INTEGER ( KIND = 4 ) :: count_,c_max,rate

    CALL SYSTEM_CLOCK(count_,rate,c_max)
    initrand=-ABS(MOD(count_,rate)+1)

  END FUNCTION INITRAND
END MODULE make_noise











 
 

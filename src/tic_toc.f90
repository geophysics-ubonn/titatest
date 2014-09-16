!!$ $Id: tic_toc.f90 1.4 2008/07/29 14:05:36 Roland Martin Exp $
MODULE tic_toc

  IMPLICIT none
  PUBLIC :: tic
  PUBLIC :: toc

  INTEGER(KIND = 4),PRIVATE         :: c2,i,l
  INTEGER(KIND = 4),PRIVATE         :: se,mi,st,ta,ms
  
CONTAINS

  SUBROUTINE tic(c1)
  ! return time in millisec into c1
    INTEGER(KIND = 4), INTENT(OUT) :: c1 ! first call -> tic
    
    CALL SYSTEM_CLOCK (c1,i)
    
  END SUBROUTINE tic
  
  
  SUBROUTINE toc(c1,csz)
  ! print elapsed time c1 to c2 in a pretty format
    INTEGER(KIND = 4), INTENT(IN) :: c1 ! first call -> tic
    CHARACTER (*),INTENT(INOUT)   :: csz
 
110 FORMAT(a,I2,':',I2.2,'.',I2.2,'h')
    
    CALL SYSTEM_CLOCK (c2,i)
    
    ms = c2-c1    ! millisec
    ms = MODULO(ms,1000)

    l = (c2-c1)/i ! sec

    ta = INT(l/24/3600) ! days
    l = l - ta*3600*24
    st = INT(l/3600) ! hours
    l = l - st*3600
    mi = INT(l/60) ! minutes
    se = l - mi*60

    WRITE (csz,110) TRIM(csz), st, mi, se
    PRINT*,TRIM(csz)

  END SUBROUTINE toc
END MODULE tic_toc

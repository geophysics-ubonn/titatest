!!$ $Id: tic_toc.f90 1.4 2008/07/29 14:05:36 Roland Martin Exp $
MODULE tic_toc

  IMPLICIT none
  PUBLIC :: tic
  PUBLIC :: toc

  INTEGER(KIND = 4),PRIVATE         :: c2,i,l
  INTEGER(KIND = 4),PRIVATE         :: se,mi,st,ta,ms
  
CONTAINS

  SUBROUTINE tic(c1)
    INTEGER(KIND = 4), INTENT(OUT) :: c1 ! first call -> tic
    
    CALL SYSTEM_CLOCK (c1,i)
    
  END SUBROUTINE tic
  
  
  SUBROUTINE toc(c1,csz)
    INTEGER(KIND = 4), INTENT(IN) :: c1 ! first call -> tic
    CHARACTER (*),INTENT(INOUT)   :: csz
 
110 FORMAT(a,I3,'d/',1X,I2,'h/',1X,I2,'m/',1X,I2,'s/',1X,I3,'ms')
    
    CALL SYSTEM_CLOCK (c2,i)
    
    ms = c2-c1    ! Gesamt Millisekunden
    ms = MODULO(ms,1000)

    l = (c2-c1)/i ! Gesamt Sekunden
    mi =INT(l/60) ! Minuten
    st =INT(mi/60) ! Stunden
    ta =INT(st/24) ! Tage
    se =l-mi*60-st*60*60-ta*60*60*24 ! Sekunden

    WRITE (csz,110)TRIM(csz),ta,st,mi,se,ms
    PRINT*,TRIM(csz)

  END SUBROUTINE toc
END MODULE tic_toc

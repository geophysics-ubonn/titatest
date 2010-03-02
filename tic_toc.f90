!!$ $Id: tic_toc.f90 1.4 2008/07/29 14:05:36 Roland Martin Exp $
MODULE tic_toc

  IMPLICIT none
  PUBLIC :: tic
  PUBLIC :: toc

  INTEGER(KIND = 4),PRIVATE,SAVE    :: c1 ! first call -> tic
  INTEGER(KIND = 4),PRIVATE         :: c2,i,l
  INTEGER(KIND = 4),PRIVATE         :: se,mi,st,ta,ms
  
CONTAINS

  SUBROUTINE tic
    
    CALL SYSTEM_CLOCK (c1,i)

  END SUBROUTINE tic


  SUBROUTINE toc

110 FORMAT(I3,'d/',1X,I2,'h/',1X,I2,'m/',1X,I2,'s/',1X,I3,'ms')

    CALL SYSTEM_CLOCK (c2,i)
    
    ms = c2-c1    ! Gesamt Millisekunden
    ms = MODULO(ms,1000)

    l = (c2-c1)/i ! Gesamt Sekunden
    mi =INT(l/60) ! Minuten
    st =INT(mi/60) ! Stunden
    ta =INT(st/24) ! Tage
    se =l-mi*60-st*60*60-ta*60*60*24 ! Sekunden

    WRITE (*,110)ta,st,mi,se,ms

  END SUBROUTINE toc
END MODULE tic_toc

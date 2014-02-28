SUBROUTINE bpar
!!!$     
!!!$     Unterprogramm zur Belegung des Parameter Vektors
!!!$     im Moment noch unspektakulaer, da elanz = manz 
!!!$     
!!!$     Copyright by Andreas Kemna
!!!$     
!!!$     Erste Version von Roland Martin                          19-Mar-2010
!!!$     
!!!$     Letzte Aenderung                                         23-Mar-2010
!!!$     
!!!$.....................................................................

  USE invmod, ONLY: par   ! fuer par
  USE sigmamod,ONLY: sigma ! fuer sigma
  USE modelmod,ONLY:manz,mnr ! fuer manz und mnr
  USE elemmod,ONLY:elanz  ! fuer elanz
  USE errmod,ONLY:fetxt,errnr ! fuer fetxt und errnr

  IMPLICIT none

!!!$.....................................................................

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsfeld
!!!$      LOGICAL,DIMENSION(:),ALLOCATABLE :: lfeld
!!!$     Indexvariablen
  INTEGER (KIND = 4) ::     i
!!!$.....................................................................


  errnr = 4

!!!$      ALLOCATE (lfeld(manz),STAT=errnr)
!!!$      
!!!$      IF (errnr/=0) THEN
!!!$         fetxt = 'Allocation problem in bpar lfeld'
!!!$         WRITE (*,'(/a/)')TRIM(fetxt)
!!!$         errnr = 97
!!!$         RETURN
!!!$      END IF
!!!$
!!!$      lfeld = .FALSE.
!!!$      PRINT*,lfeld
!!!$     Parametervektor belegen
!!!$      do i=1,elanz
!!!$         if (.not.lfeld(j)) then
!!!$            lfeld(j) = .true.
!!!$            par(j)   = LOG(sigma(i))
!!!$         end if
!!!$      end do

  do i=1,elanz
     par(mnr(i))   = LOG(sigma(i))
  end do

  errnr = 0

end SUBROUTINE bpar

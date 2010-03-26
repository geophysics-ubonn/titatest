      SUBROUTINE bpar
c     
c     Unterprogramm zur Beleguing des Paramter vektors
c     im Moment noch unspektakulaer, da elanz = manz 
c     
c     Copyright by Andreas Kemna
c     
c     Erste Version von Roland Martin                          19-Mar-2010
c     
c     Letzte Aenderung                                         23-Mar-2010
c     
c.....................................................................
      IMPLICIT none
      INCLUDE 'parmax.fin' ! fuer maximale zahlenwerte..
      INCLUDE 'sigma.fin' ! fuer sigma
      INCLUDE 'model.fin' ! fuer manz und mnr
      INCLUDE 'elem.fin'  ! fuer elanz
      INCLUDE 'inv.fin'   ! fuer par
      INCLUDE 'err.fin' ! fuer fetxt und errnr
c.....................................................................

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Hilfsfeld
c      LOGICAL,DIMENSION(:),ALLOCATABLE :: lfeld
c     Indexvariablen
      integer         * 4     i
c.....................................................................


      errnr = 4

c$$$      ALLOCATE (lfeld(manz),STAT=errnr)
c$$$      
c$$$      IF (errnr/=0) THEN
c$$$         fetxt = 'Allocation problem in bpar lfeld'
c$$$         WRITE (*,'(/a/)')TRIM(fetxt)
c$$$         errnr = 97
c$$$         RETURN
c$$$      END IF
c$$$
c$$$      lfeld = .FALSE.
c$$$      PRINT*,lfeld
c     Parametervektor belegen
c$$$      do i=1,elanz
c$$$         if (.not.lfeld(j)) then
c$$$            lfeld(j) = .true.
c$$$            par(j)   = CDLOG(sigma(i))
c$$$         end if
c$$$      end do

      do i=1,elanz
         par(mnr(i))   = CDLOG(sigma(i))
      end do

      errnr = 0

      end

      SUBROUTINE bsigma
c     
c     Unterprogramm zur Belegung des Modell vektors (sigma)
c     mit dem verbesserten Modell aus par
c     im Moment noch unspektakulaer, da elanz = manz 
c     
c     Copyright by Andreas Kemna
c     
c     Erste Version von Roland Martin                          19-Mar-2010
c     
c     Letzte Aenderung                                         23-Mar-2010
c     
c.....................................................................

      USE invmod   ! fuer par

      IMPLICIT none

      INCLUDE 'parmax.fin' ! fuer maximale zahlenwerte..
      INCLUDE 'sigma.fin' ! fuer sigma
      INCLUDE 'model.fin' ! fuer manz und mnr
      INCLUDE 'elem.fin'  ! fuer elanz
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

      do i=1,elanz
         sigma(i) = CDEXP(par(mnr(i)))
      end do

      errnr = 0

      end

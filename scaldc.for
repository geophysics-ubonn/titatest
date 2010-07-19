      subroutine scaldc()
      
c     Unterprogramm skaliert 'adc' und 'bdc' und liefert die Skalierungs-
c     faktoren im Vektor 'fak'.

c     ( Vgl. Subroutine 'SCALBNDN' in Schwarz (1991) )

c     Andreas Kemna                                            11-Oct-1993
c     Letzte Aenderung   07-Mar-2003

c.....................................................................

      USE alloci
      USE femmod

      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'
      INCLUDE 'elem.fin'

c.....................................................................

c     Hilfsvariablen
      integer         * 4     idi,i0
      integer         * 4     ja
      real            * 8     dum

c     Indexvariablen
      integer         * 4     i,j

c.....................................................................

      do i=1,sanz

         idi = i*(mb+1)
         dum = adc(idi)

         if (dum.le.0d0) then
            WRITE (fetxt,'(a,I6,A,I6)')
     1           'scaldc',i,'idi',idi
            errnr = 27
            RETURN
         end if

         adc(idi) = 1d0
         fak(i)   = 1d0 / dsqrt(dum)
         bdc(i)   = bdc(i) * fak(i)

         if (i.eq.1) CYCLE

         i0 = i*mb
         ja = max0(1,i-mb)

         do j=ja,i-1
            adc(i0+j) = adc(i0+j)*fak(i)*fak(j)
         end do

      END DO

      errnr = 0

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

      end

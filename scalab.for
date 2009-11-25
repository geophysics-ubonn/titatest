      subroutine scalab()
      
c     Unterprogramm skaliert 'a' und 'b' und liefert die Skalierungs-
c     faktoren im Vektor 'fak'.

c     ( Vgl. Subroutine 'SCALBNDN' in Schwarz (1991) )

c     Andreas Kemna                                            11-Oct-1993
c     Letzte Aenderung   07-Mar-2003

c.....................................................................

      USE alloci
      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'fem.fin'

c.....................................................................

c     Hilfsvariablen
      integer         * 4     idi,i0
      integer         * 4     ja
      real            * 8     dum

c     Indexvariablen
      integer         * 4     i,j

c.....................................................................

      do 30 i=1,sanz

         idi = i*(mb+1)
         dum = cdabs(a(idi))

         if (dum.le.0d0) then
            WRITE (fetxt,*)
     1           'scalab idi',idi,'i',i
            errnr = 27
            goto 1000
         end if

         a(idi) = a(idi) / dcmplx(dum)
         fak(i) = 1d0 / dsqrt(dum)
         b(i)   = b(i) * dcmplx(fak(i))

         if (i.eq.1) goto 30

         i0 = i*mb
         ja = max0(1,i-mb)

         do 20 j=ja,i-1
            a(i0+j) = a(i0+j) * dcmplx(fak(i)*fak(j))
 20      continue

 30   continue

      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 1000 return

      end

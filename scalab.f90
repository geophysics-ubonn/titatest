      subroutine scalab()
      
!!!$     Unterprogramm skaliert 'a' und 'b' und liefert die Skalierungs-
!!!$     faktoren im Vektor 'fak'.

!!!$     ( Vgl. Subroutine 'SCALBNDN' in Schwarz (1991) )

!!!$     Andreas Kemna                                            11-Oct-1993
!!!$     Letzte Aenderung   07-Mar-2003

!!!$.....................................................................

      USE alloci
      USE femmod
      USE elemmod
      USE errmod

      IMPLICIT none


!!!$.....................................................................

!!!$     Hilfsvariablen
      integer         * 4     idi,i0
      integer         * 4     ja
      real            * 8     dum

!!!$     Indexvariablen
      integer         * 4     i,j

!!!$.....................................................................

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

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

 1000 return

      end

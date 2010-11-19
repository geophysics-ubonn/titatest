      subroutine choldc()
      
c     Cholesky-Zerlegung der positiv definiten Matrix 'adc'; erfolgt auf dem
c     Platz von 'adc', d.h. bei Auftreten eines Fehlers ist gegebene Matrix
c     'adc' zerstoert.

c     ( Vgl. Subroutine 'CHOBNDN' in Schwarz (1991) )

c     Andreas Kemna                                            11-Oct-1993
c     Letzte Aenderung   07-Mar-2003

c.....................................................................

      USE alloci
      USE elemmod
      USE errmod

      IMPLICIT none


c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Hilfsvariablen
      INTEGER (KIND = 4)  ::     idi,i0,ij,j0
      INTEGER (KIND = 4)  ::     m1,fi
      REAL (KIND(0D0))    ::     s

c     Indexvariablen
      INTEGER (KIND = 4)  ::     i,j,k

c.....................................................................

      m1 = mb+1

      do 30 i=1,sanz

         idi = i*m1
         fi  = max0(1,i-mb)
         i0  = idi-i

         do 20 j=fi,i

            ij = i0+j
            j0 = j*mb
            s  = adc(ij)

            do 10 k=fi,j-1
               s = s - adc(i0+k)*adc(j0+k)
 10         continue

            if (j.lt.i) then

               adc(ij) = s / adc(j*m1)

            else

               if (s.le.0d0) then
                  fetxt = ' '
                  errnr = 28
                  goto 1000
               end if

               adc(idi) = dsqrt(s)

            end if

 20      continue

 30   continue

      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 1000 return

      end

      subroutine randb2()

c     Unterprogramm modifiziert die Matrix 'a' (Bandbreite 'mb') und den
c     Konstantenvektor 'b' zur Beruecksichtigung der Dirichletschen Rand-
c     bedingungen ('rwdanz' Randwerte 'rwd(rwdmax)' mit zugeh. Knotennummern
c     'rwdnr(smax)').

c     ( Vgl. Subroutine 'RBSTBNDN' in Schwarz (1991) )

c     Andreas Kemna                                            12-Feb-1993
c     Letzte Aenderung   15-Jul-2007

c.....................................................................

      USE alloci
      USE femmod

      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'randb.fin'
      INCLUDE 'elem.fin'

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Hilfsvariablen
      integer         * 4     m1
      real            * 8     rwertdc
      complex         * 16    rwert

c     Indexvariablen
      integer         * 4     ir,i,j,k,
     1     idk,ia,ki,ji

c.....................................................................

      if (rwdanz.eq.0) return

      m1 = mb + 1

      do 50 ir=1,rwdanz
         k      = rwdnr(ir)

         if (ldc) then
            rwertdc  = rwddc(ir)
            bdc(k)   = -rwertdc
            idk      = k*m1
            adc(idk) = 1d0
         else
            rwert  = rwd(ir)
            b(k)   = -rwert
            idk    = k*m1
            a(idk) = dcmplx(1d0)
         end if

         if (k.eq.1) goto 30

         ia = max0(1,mb+2-k)

         do 20 i=ia,mb
            j  = k+i-m1
            ki = idk+i-m1

            if (ldc) then
               bdc(j)  = bdc(j) + rwertdc * adc(ki)
               adc(ki) = (0d0)
            else
               b(j)  = b(j) + rwert * a(ki)
               a(ki) = dcmplx(0d0)
            end if
 20      continue

 30      if (k.eq.sanz) goto 50

         ia = max0(1,k-sanz+m1)

         do 40 i=ia,mb
            j  = k-i+m1
            ji = (j-1)*m1+i

            if (ldc) then
               bdc(j)  = bdc(j) + rwertdc * adc(ji)
               adc(ji) = (0d0)
            else
               b(j)  = b(j) + rwert * a(ji)
               a(ji) = dcmplx(0d0)
            end if
 40      continue

 50   continue

      return
      end

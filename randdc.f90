      subroutine randdc()
      
!!!$     Unterprogramm modifiziert die Matrix 'adc' (Bandbreite 'mb') und den
!!!$     Konstantenvektor 'bdc' gemaess homogener Dirichletscher Randbedingungen.

!!!$     ( Vgl. Subroutine 'RBSTBNDN' in Schwarz (1991) )

!!!$     Andreas Kemna                                            12-Feb-1993
!!!$     Letzte Aenderung   13-Nov-1997

!!!$.....................................................................

      USE alloci
      USE femmod
      USE elemmod

      IMPLICIT none

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Aktuelle Elementnummer
      integer         * 4     iel

!!!$     Aktueller Elementtyp
      integer         * 4     ntyp

!!!$     Anzahl der Knoten im aktuellen Elementtyp
      integer         * 4     nkel

!!!$     Index-/Hilfsvariablen
      integer         * 4     m1,ir,i,j,k,ia,ki,i2,j2,idk,ji
      
!!!$.....................................................................

      m1  = mb+1
      iel = 0

      do i2=1,typanz
         ntyp = typ(i2)
         nkel = selanz(i2)

         do 30 j2=1,nelanz(i2)
            iel = iel+1

            if (ntyp.ne.13) goto 30

            do ir=1,nkel
               k      = nrel(iel,ir)
               bdc(k) = 0d0

               idk      = k*m1
               adc(idk) = 1d0

               if (k.eq.1) goto 10

               ia = max0(1,mb+2-k)

               do i=ia,mb
                  ki      = idk+i-m1
                  adc(ki) = 0d0
               end do

 10            if (k.eq.sanz) goto 20

               ia = max0(1,k-sanz+m1)

               do i=ia,mb
                  j       = k-i+m1
                  ji      = (j-1)*m1+i
                  adc(ji) = 0d0
               end do

 20            continue
            end do
 30      continue
      end do

      return
      end

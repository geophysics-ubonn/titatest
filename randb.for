      subroutine randb()
      
c     Unterprogramm modifiziert die Matrix 'a' (Bandbreite 'mb') und den
c     Konstantenvektor 'b' gemaess homogener Dirichletscher Randbedingungen.

c     ( Vgl. Subroutine 'RBSTBNDN' in Schwarz (1991) )

c     Andreas Kemna                                            12-Feb-1993
c     Letzte Aenderung   05-Nov-1997

c.....................................................................

      USE alloci
      USE femmod

      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'elem.fin'

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Aktuelle Elementnummer
      integer         * 4     iel

c     Aktueller Elementtyp
      integer         * 4     ntyp

c     Anzahl der Knoten im aktuellen Elementtyp
      integer         * 4     nkel

c     Index-/Hilfsvariablen
      integer         * 4     m1,ir,i,j,k,ia,ki,i2,j2,idk,ji
      
c.....................................................................

      m1  = mb+1
      iel = 0

      do i2=1,typanz
         ntyp = typ(i2)
         nkel = selanz(i2)

         do 30 j2=1,nelanz(i2)
            iel = iel+1

            if (ntyp.ne.13) goto 30

            do ir=1,nkel
               k    = nrel(iel,ir)
               b(k) = dcmplx(0d0)

               idk    = k*m1
               a(idk) = dcmplx(1d0)

               if (k.eq.1) goto 10

               ia = max0(1,mb+2-k)

               do i=ia,mb
                  ki    = idk+i-m1
                  a(ki) = dcmplx(0d0)
               end do

 10            if (k.eq.sanz) goto 20

               ia = max0(1,k-sanz+m1)

               do i=ia,mb
                  j     = k-i+m1
                  ji    = (j-1)*m1+i
                  a(ji) = dcmplx(0d0)
               end do

 20            continue
            end do
 30      continue
      end do

      return
      end

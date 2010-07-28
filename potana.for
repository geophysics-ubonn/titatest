      subroutine potana(l,k)

c     Unterprogramm zur Berechnung der analytischen Loesung der
c     Helmholtzgleichung fuer einen homogenen Halbraum (fuer Einheitsstrom !).

c     Andreas Kemna                                            04-Jan-1996
c     Letzte Aenderung   11-Nov-1997
      
c.....................................................................

      USE femmod
      USE sigmamod
      USE electrmod
      USE elemmod
      USE wavenmod

      IMPLICIT none

c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Elektrodenindex
      integer         * 4     l

c     Wellenzahlindex
      integer         * 4     k

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Hilfsfunction
      real            * 8     bessk0

c     Hilfsvariablen
      real            * 8     xk1,yk1,
     1     xk2,yk2,
     1     x21,y21m,y21p,
     1     rm,rp,potmax,dum
      complex         * 16    dum2
      integer         * 4     idum,j

c     Pi
      real            * 8     pi

c.....................................................................

      pi = dacos(-1d0)

      xk1    = sx(snr(enr(l)))
      yk1    = sy(snr(enr(l)))
      idum   = enr(l)
      potmax = 0d0

      do j=1,idum-1
         xk2  = sx(snr(j))
         yk2  = sy(snr(j))
         x21  = xk2-xk1
         y21m = yk2-yk1
         y21p = yk2+yk1
         rm   = dsqrt(x21*x21+y21m*y21m)
         rp   = dsqrt(x21*x21+y21p*y21p)

         dum     = bessk0(rm*kwn(k)) + bessk0(rp*kwn(k))
         potmax  = dmax1(potmax,dum)
         pota(j) = dcmplx(dum)
      end do

      do j=idum+1,sanz
         xk2  = sx(snr(j))
         yk2  = sy(snr(j))
         x21  = xk2-xk1
         y21m = yk2-yk1
         y21p = yk2+yk1
         rm   = dsqrt(x21*x21+y21m*y21m)
         rp   = dsqrt(x21*x21+y21p*y21p)

         dum     = bessk0(rm*kwn(k)) + bessk0(rp*kwn(k))
         potmax  = dmax1(potmax,dum)
         pota(j) = dcmplx(dum)
      end do

c     Endlichen Wert fuer Singularitaet vorgeben (beeinflusst nur
c     berechnete Potentialwerte in direkter Umgebung des Stromknotens !)
c     ak
      pota(idum) = dcmplx(5d0*potmax)

c     Potentialwerte skalieren (fuer Einheitsstrom !)
      dum2 = dcmplx(5d-1/pi) / sigma0

      do j=1,sanz
         pota(j) = pota(j) * dum2
      end do

      return
      end

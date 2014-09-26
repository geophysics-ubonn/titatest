subroutine potana(l,k,my_pota)

!!!$     Unterprogramm zur Berechnung der analytischen Loesung der
!!!$     Helmholtzgleichung fuer einen homogenen Halbraum (fuer Einheitsstrom !).

!!!$     Andreas Kemna                                            04-Jan-1996
!!!$     Letzte Aenderung   11-Nov-1997

!!!$.....................................................................
  use alloci, only: prec
  USE sigmamod
  USE electrmod
  USE elemmod, ONLY : sanz, snr, sx, sy
  USE wavenmod

  IMPLICIT none

!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Elektrodenindex
  INTEGER (KIND=4)  :: l

!!!$     Wellenzahlindex
  INTEGER (KIND=4)  :: k

!!$ Analytische berechnete Potentialwerte 
  COMPLEX (prec), DIMENSION(sanz) :: my_pota

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsfunction
  REAL (prec)  :: bessk0

!!!$     Hilfsvariablen
  REAL (prec)  :: xk1,yk1,xk2,yk2,x21,y21m,y21p,rm,rp,potmax,dum
  COMPLEX (prec)  ::    dum2
  INTEGER (KIND=4)  ::  idum,j

!!!$     Pi
  REAL (prec)  ::  pi

!!!$.....................................................................

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
     rm   = SQRT(x21*x21+y21m*y21m)
     rp   = SQRT(x21*x21+y21p*y21p)

     dum     = bessk0(rm*kwn(k)) + bessk0(rp*kwn(k))
     potmax  = (MAX(potmax,dum))
     my_pota(j) = dCMPLX(dum)
  end do

  do j=idum+1,sanz
     xk2  = sx(snr(j))
     yk2  = sy(snr(j))
     x21  = xk2-xk1
     y21m = yk2-yk1
     y21p = yk2+yk1
     rm   = SQRT(x21*x21+y21m*y21m)
     rp   = SQRT(x21*x21+y21p*y21p)

     dum     = bessk0(rm*kwn(k)) + bessk0(rp*kwn(k))
     potmax  = (MAX(potmax,dum))
     my_pota(j) = dCMPLX(dum)
  end do

!!!$     Endlichen Wert fuer Singularitaet vorgeben (beeinflusst nur
!!!$     berechnete Potentialwerte in direkter Umgebung des Stromknotens !)
!!!$     ak
  my_pota(idum) = dCMPLX(5d0*potmax)

!!!$     Potentialwerte skalieren (fuer Einheitsstrom !)
  dum2 = dCMPLX(5d-1/pi) / sigma0

  my_pota = my_pota * dum2
end subroutine potana

function beta(nelec,k)

!!!$     Function zur Berechnung der 'mixed boundary conditions'.

!!!$     Andreas Kemna                                         20-Dec-1993
!!!$     Letzte Aenderung                                      20-Nov-2009

!!!$.....................................................................
use alloci, only: prec
  USE electrmod
  USE elemmod
  USE wavenmod
  USE errmod

  IMPLICIT none


!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Aktuelle Elektrodennummer
  INTEGER (KIND = 4)  ::     nelec

!!!$     Aktueller Wellenzahlindex
  INTEGER (KIND = 4)  ::     k

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:
real(prec) beta
!!!$     Abstand Randelement - Quelle/Spiegelquelle
  REAL (prec)    ::     r1m,r1p

!!!$     Hilfsfunctions
  REAL (prec)    ::     bessk0,bessk1

!!!$     Hilfsvariablen
  REAL (prec)    ::     cthetm,cthetp
  REAL (prec)    ::     xs,ys,xr,yr,x3,y3m,y3p,&
       x4,y4,r2,bk0m,bk0p,bk1m,bk1p

!!!$.....................................................................

!!!$     Koordinaten der Stromelektrode bestimmen
  xs = sx(snr(enr(nelec)))
  ys = sy(snr(enr(nelec)))

!!!$     Koordinaten des aktuellen Randelements bestimmen
  xr = (xk(1)+xk(2)) / 2d0
  yr = (yk(1)+yk(2)) / 2d0

!!!$     Abstand Randelement - Quelle/Spiegelquelle bestimmen
  x3  = xr - xs
  y3m = (yr - sytop) - (ys - sytop)
  y3p = (yr - sytop) + (ys - sytop)

  r1m = SQRT(x3*x3 + y3m*y3m)
  r1p = SQRT(x3*x3 + y3p*y3p)

!!!$     Hilfsvektor bestimmen
  x4 = yk(1) - yk(2)
  y4 = xk(2) - xk(1)

  r2 = SQRT(x4*x4 + y4*y4)

!!!$     Ggf. Fehlermeldung
  if (r1m.lt.1d-12.or.r1p.lt.1d-12.or.r2.lt.1d-12) then
     fetxt = ' '
     errnr = 35
     goto 1000
  end if

!!!$     Kosinus von Theta bestimmen
  cthetm = (x3*x4 + y3m*y4) / r1m / r2
  cthetp = (x3*x4 + y3p*y4) / r1p / r2

!!!$     'beta' bestimmen
  bk0m = bessk0(r1m*kwn(k))
  bk0p = bessk0(r1p*kwn(k))

  bk1m = bessk1(r1m*kwn(k))
  bk1p = bessk1(r1p*kwn(k))

  if (bk0m+bk0p.eq.0d0) then
     beta = (cthetm+cthetp) / 2d0
  else
     beta = (cthetm*bk1m+cthetp*bk1p) / (bk0m+bk0p)
  end if

  beta = kwn(k) * beta
  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

1000 return

end function beta

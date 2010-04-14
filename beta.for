      real * 8 function beta(nelec,k,xk,yk)

c     Function zur Berechnung der 'mixed boundary conditions'.

c     Andreas Kemna                                         20-Dec-1993
c     Letzte Aenderung                                      20-Nov-2009

c.....................................................................
      IMPLICIT none
      INCLUDE 'parmax.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'electr.fin'
      INCLUDE 'waven.fin'
      INCLUDE 'err.fin'

c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Aktuelle Elektrodennummer
      integer         * 4     nelec

c     Aktueller Wellenzahlindex
      integer         * 4     k

c     x-Koordinaten der Eckknotenpunkte
      real            * 8     xk(selmax)

c     y-Koordinaten der Eckknotenpunkte
      real            * 8     yk(selmax)

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Abstand Randelement - Quelle/Spiegelquelle
      real            * 8     r1m,r1p

c     Hilfsfunctions
      real            * 8     bessk0,
     1     bessk1

c     Hilfsvariablen
      real            * 8     cthetm,cthetp
      real            * 8     xs,ys,xr,yr,
     1     x3,y3m,y3p,
     1     x4,y4,r2,
     1     bk0m,bk0p,
     1     bk1m,bk1p

c.....................................................................

c     Koordinaten der Stromelektrode bestimmen
      xs = sx(snr(enr(nelec)))
      ys = sy(snr(enr(nelec)))

c     Koordinaten des aktuellen Randelements bestimmen
      xr = (xk(1)+xk(2)) / 2d0
      yr = (yk(1)+yk(2)) / 2d0

c     Abstand Randelement - Quelle/Spiegelquelle bestimmen
      x3  = xr - xs
      y3m = (yr - sytop) - (ys - sytop)
      y3p = (yr - sytop) + (ys - sytop)

      r1m = dsqrt(x3*x3 + y3m*y3m)
      r1p = dsqrt(x3*x3 + y3p*y3p)

c     Hilfsvektor bestimmen
      x4 = yk(1) - yk(2)
      y4 = xk(2) - xk(1)

      r2 = dsqrt(x4*x4 + y4*y4)

c     Ggf. Fehlermeldung
      if (r1m.lt.1d-12.or.r1p.lt.1d-12.or.r2.lt.1d-12) then
         fetxt = ' '
         errnr = 35
         goto 1000
      end if

c     Kosinus von Theta bestimmen
      cthetm = (x3*x4 + y3m*y4) / r1m / r2
      cthetp = (x3*x4 + y3p*y4) / r1p / r2

c     'beta' bestimmen
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

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 1000 return

      end

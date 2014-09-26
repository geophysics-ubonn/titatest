subroutine refsig()

!!!$     Unterprogramm zum Bestimmen der Referenzleitfaehigkeit.

!!!$     Andreas Kemna                                            29-Feb-1996
!!!$     Letzte Aenderung   08-Nov-1997

!!!$.....................................................................
use alloci, only: prec
  USE sigmamod
  USE elemmod,only:sx,sy,snr,typanz,nelanz,typ,nrel

  IMPLICIT none

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Indexvariablen
  INTEGER (KIND = 4 ) ::     i,j

!!!$     Hilfsvariablen
  INTEGER (KIND = 4) ::      iel
  REAL (prec)   ::     xk,yk,ax,ay,bx,by,cx,cy,dum,area

!!!$.....................................................................

  iel    = 0
  area   = 0d0
  sigma0 = dCMPLX(0d0)

  do i=1,typanz
     do j=1,nelanz(i)
        iel = iel + 1

        if (typ(i).gt.10) CYCLE

        xk = sx(snr(nrel(iel,1)))
        yk = sy(snr(nrel(iel,1)))
        ax = sx(snr(nrel(iel,2)))
        ay = sy(snr(nrel(iel,2)))
        bx = sx(snr(nrel(iel,3)))
        by = sy(snr(nrel(iel,3)))
        ax = ax - xk
        ay = ay - yk
        bx = bx - xk
        by = by - yk

        if (typ(i).eq.3) then
           dum = ABS(ax*by-ay*bx) / 2d0
        else if (typ(i).eq.5.or.typ(i).eq.8) then
           cx  = sx(snr(nrel(iel,4)))
           cy  = sy(snr(nrel(iel,4)))
           cx  = cx - xk
           cy  = cy - yk
           dum = (ABS(ax*by-ay*bx) + ABS(bx*cy-by*cx)) / 2d0
        end if

        sigma0 = sigma0 + dCMPLX(dum)*LOG(sigma(iel))
        area   = area   + dum

     end do
  end do

  sigma0 = EXP(sigma0/dCMPLX(area))

  return
end subroutine refsig

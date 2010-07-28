      subroutine refsig()
      
c     Unterprogramm zum Bestimmen der Referenzleitfaehigkeit.

c     Andreas Kemna                                            29-Feb-1996
c     Letzte Aenderung   08-Nov-1997

c.....................................................................

      USE sigmamod
      USE elemmod

      IMPLICIT none

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Indexvariablen
      integer         * 4     i,j

c     Hilfsvariablen
      integer         * 4     iel
      real            * 8     xk,yk,
     1     ax,ay,bx,by,cx,cy,
     1     dum,area

c.....................................................................

      iel    = 0
      area   = 0d0
      sigma0 = dcmplx(0d0)

      do i=1,typanz
         do j=1,nelanz(i)
            iel = iel + 1

            if (typ(i).gt.10) goto 10

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
               dum = dabs(ax*by-ay*bx) / 2d0
            else if (typ(i).eq.5.or.typ(i).eq.8) then
               cx  = sx(snr(nrel(iel,4)))
               cy  = sy(snr(nrel(iel,4)))
               cx  = cx - xk
               cy  = cy - yk
               dum = (dabs(ax*by-ay*bx) + dabs(bx*cy-by*cx)) / 2d0
            end if

            sigma0 = sigma0 + dcmplx(dum)*cdlog(sigma(iel))
            area   = area   + dum

 10         continue
         end do
      end do

      sigma0 = cdexp(sigma0/dcmplx(area))

      return
      end

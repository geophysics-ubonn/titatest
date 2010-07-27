      subroutine cjggra()

c     Unterprogramm berechnet Modellverbesserung mittels konjugierter
c     Gradienten.

c     Andreas Kemna                                            01-Mar-1996
c     Letzte Aenderung   04-Jul-2009
c.....................................................................

      USE invmod
      USE cjgmod

      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'model.fin'
      INCLUDE 'konv.fin'
      INCLUDE 'err.fin'

c.....................................................................
c     PROGRAMMINTERNE PARAMETER:

c     Skalare
      complex         * 16    beta
      real            * 8     alpha,
     1     dr,dr0,dr1

c     Hilfsvariablen
      integer         * 4     j,k

c.....................................................................
      
      ALLOCATE (rvec(manz),pvec(manz),ap(manz),
     1     stat=errnr)
      IF (errnr /= 0) THEN
         fetxt = 'Error memory allocation rvec in cjggdc'
         errnr = 94
         RETURN
      END IF

      dpar = 0.
      rvec = bvec
      pvec = 0.

      do k=1,ncgmax
         ncg = k-1

         dr = 0d0
         do j=1,manz
            dr = dr + dble(dconjg(rvec(j))*rvec(j))
         end do

         if (k.eq.1) then
            dr0  = dr*eps
            beta = dcmplx(0d0)
         else
            if (dr.le.dr0) goto 10

c     Fletcher-Reeves-Version
            beta = dcmplx(dr/dr1)

c     akc Polak-Ribiere-Version
c     ak                beta = dcmplx(0d0)
c     ak                do j=1,manz
c     ak                    beta = beta + dconjg(bvec(j))*rvec(j)
c     ak                end do
c     ak                beta = beta*dcmplx(-alpha/dr1)
         end if

         do j=1,manz
            pvec(j) = rvec(j) + beta*pvec(j)
         end do
         
         IF (ltri == 0) THEN
            CALL bp
         ELSE IF (ltri == 1.OR.ltri == 2.OR.
     1           (ltri > 4 .AND. ltri < 15)) THEN
            call bptri
         ELSE IF (ltri == 3.OR.ltri == 4) THEN
            CALL bplma
         ELSE IF (ltri == 15) THEN
            call bpsto
         END IF

         dr1 = 0d0
         do j=1,manz
            dr1 = dr1 + dble(dconjg(pvec(j))*bvec(j))
         end do

         alpha = dr/dr1

         do j=1,manz
            dpar(j) = dpar(j) + dcmplx(alpha)*pvec(j)
            rvec(j) = rvec(j) - dcmplx(alpha)*bvec(j)
         end do

         dr1 = dr

c     Residuum speichern
         cgres(k+1) = real(eps*dr/dr0)
      end do

      ncg = ncgmax

c     Anzahl an CG-steps speichern
 10   cgres(1) = real(ncg)

      DEALLOCATE (rvec,pvec,ap)

      end

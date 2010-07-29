      subroutine cjggra()

c     Unterprogramm berechnet Modellverbesserung mittels konjugierter
c     Gradienten.

c     Andreas Kemna                                            01-Mar-1996
c     Letzte Aenderung   04-Jul-2009
c.....................................................................

      USE invmod
      USE cjgmod
      USE modelmod, ONLY : manz
      USE datmod, ONLY : nanz

      IMPLICIT none

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
      
      ALLOCATE (rvec(manz),pvec(manz),ap(nanz),
     1     stat=errnr)
      IF (errnr /= 0) THEN
         fetxt = 'Error memory allocation rvec in cjggdc'
         errnr = 94
         RETURN
      END IF

      dpar = 0.
      rvec = bvec
      pvec = 0.

      fetxt = 'CG iteration'

      do k=1,ncgmax

         ncg = k-1

c$$$         dr = 0d0
c$$$         do j=1,manz
c$$$            dr = dr + dble(dconjg(rvec(j))*rvec(j))
c$$$         end do

         dr = DOT_PRODUCT(DCONJG(rvec),rvec)

         if (k.eq.1) then
            dr0  = dr*eps
            beta = dcmplx(0d0)
         else
c     Fletcher-Reeves-Version
            beta = dcmplx(dr/dr1)
c     akc Polak-Ribiere-Version
c     ak                beta = dcmplx(0d0)
c     ak                do j=1,manz
c     ak                    beta = beta + dconjg(bvec(j))*rvec(j)
c     ak                end do
c     ak                beta = beta*dcmplx(-alpha/dr1)
         END IF

         WRITE (*,'(a,t40,I5,t55,G10.4,t70,G10.4)',ADVANCE='no')
     1        ACHAR(13)//TRIM(fetxt),k,dr,dr0
         
         if (dr.le.dr0) goto 10

         pvec = rvec + beta * pvec

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

         dr1 = DOT_PRODUCT(DCONJG(pvec),bvec)
c$$$
c$$$         dr1 = 0d0
c$$$         do j=1,manz
c$$$            dr1 = dr1 + dble(dconjg(pvec(j))*bvec(j))
c$$$         end do

         alpha = dr/dr1

         dpar = dpar + DCMPLX(alpha) * pvec
         rvec = rvec - DCMPLX(alpha) * bvec

         dr1 = dr

c     Residuum speichern
         cgres(k+1) = real(eps*dr/dr0)
      end do

      ncg = ncgmax

c     Anzahl an CG-steps speichern
 10   cgres(1) = real(ncg)

      DEALLOCATE (rvec,pvec,ap)

      end

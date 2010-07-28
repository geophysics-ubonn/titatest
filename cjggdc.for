      subroutine cjggdc()

c     Unterprogramm berechnet Modellverbesserung mittels konjugierter
c     Gradienten.

c     Andreas Kemna                                            01-Mar-1996
c     Letzte Aenderung                                         29-Jul-2009
c.....................................................................

      USE invmod
      USE cjgmod
      USE modelmod

      IMPLICIT none

      INCLUDE 'konv.fin'
      INCLUDE 'err.fin'

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Skalare
      real            * 8     beta,alpha,
     1     dr,dr0,dr1

c     Hilfsvariablen
      integer         * 4     j,k

c.....................................................................
      ALLOCATE (rvecdc(manz),pvecdc(manz),apdc(manz),
     1     bvecdc(manz),stat=errnr)
      IF (errnr /= 0) THEN
         fetxt = 'Error memory allocation rvec in cjggdc'
         errnr = 94
         RETURN
      END IF
       
      if (lip) then
         bvecdc = dimag(bvec)
      else
         bvecdc = dble(bvec)
      end if
      dpar = DCMPLX(0D0)
      rvecdc = bvecdc
      pvecdc = 0.

      do k=1,ncgmax
         ncg = k-1

         dr = DOT_PRODUCT(rvecdc,rvecdc)

         if (k.eq.1) then
            dr0  = dr*eps
            beta = 0d0
         else
            beta = dr/dr1
         end if

         if (dr.le.dr0) goto 10

         pvecdc = rvecdc + beta * pvecdc

         IF (ltri == 0) THEN
            CALL bpdc

         ELSE IF (ltri == 1.OR.ltri == 2.OR.
     1           (ltri > 4 .AND. ltri < 15)) THEN
            CALL bpdctri

         ELSE IF (ltri == 3.OR.ltri == 4) THEN
            CALL bpdclma

         ELSE IF (ltri == 15) THEN
            CALL bpdcsto

         END IF

         dr1 = DOT_PRODUCT(pvecdc,bvecdc)

         alpha = dr/dr1

         do j=1,manz
            dpar(j) = dpar(j) + dcmplx(alpha*pvecdc(j))
            rvecdc(j) = rvecdc(j) - alpha*bvecdc(j)
         end do

         dr1 = dr

c     Residuum speichern
         cgres(k+1) = real(eps*dr/dr0)
      end do

      ncg = ncgmax

c     Anzahl an CG-steps speichern
 10   cgres(1) = real(ncg)

      DEALLOCATE (pvecdc,rvecdc,apdc,bvecdc)

      end

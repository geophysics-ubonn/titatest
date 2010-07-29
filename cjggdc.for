      subroutine cjggdc()

c     Unterprogramm berechnet Modellverbesserung mittels konjugierter
c     Gradienten.

c     Andreas Kemna                                            01-Mar-1996
c     Letzte Aenderung                                         29-Jul-2009
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
      REAL(KIND(0D0)) :: beta,alpha,dr,dr0,dr1

c     Hilfsvariablen
      integer         * 4     j,k

c.....................................................................
      ALLOCATE (rvecdc(manz),pvecdc(manz),apdc(nanz),
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

      dpar = DCMPLX(0.)
      rvecdc = bvecdc
      pvecdc = 0.

      fetxt = 'CG iteration'

      do k=1,ncgmax

         ncg = k-1

         dr = DOT_PRODUCT(rvecdc,rvecdc)

         if (k.eq.1) then
            dr0  = dr*eps
            beta = 0d0
         else
            beta = dr/dr1
         end if

         WRITE (*,'(a,t40,I5,t55,G10.4,t70,G10.4)',ADVANCE='no')
     1        ACHAR(13)//TRIM(fetxt),k,dr,dr0

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

         dpar = dpar + DCMPLX(alpha) * DCMPLX(pvecdc)
         rvecdc = rvecdc - alpha * bvecdc

c rm update speichern
         dr1 = dr

c     Residuum speichern
         cgres(k+1) = real(eps*dr/dr0)
      end do

      ncg = ncgmax

c     Anzahl an CG-steps speichern
 10   cgres(1) = real(ncg)

      DEALLOCATE (pvecdc,rvecdc,apdc,bvecdc)

      end

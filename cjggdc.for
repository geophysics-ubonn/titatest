        subroutine cjggdc(bvec)

c Unterprogramm berechnet Modellverbesserung mittels konjugierter
c Gradienten.

c Andreas Kemna                                            01-Mar-1996
c                                       Letzte Aenderung   09-Jan-1998
         
c.....................................................................

        INCLUDE 'parmax.fin'
        INCLUDE 'model.fin'
        INCLUDE 'inv.fin'
        INCLUDE 'konv.fin'

c.....................................................................

c EIN-/AUSGABEPARAMETER:

c Konstantenvektor
        complex         * 16    bvec(mmax)

c.....................................................................

c PROGRAMMINTERNE PARAMETER:

c Vektoren
        real            * 8     rvec(mmax),pvec(mmax)
        real            * 8     bvecdc(mmax)

c Skalare
        real            * 8     beta,alpha,
     1                          dr,dr0,dr1

c Hilfsvariablen
        integer         * 4     j,k

c.....................................................................

        do j=1,manz
            if (lip) then
                bvecdc(j) = dimag(bvec(j))
            else
                bvecdc(j) = dble(bvec(j))
            end if

            dpar(j) = dcmplx(0d0)
            rvec(j) = bvecdc(j)
            pvec(j) = 0d0
        end do

        do k=1,ncgmax
            ncg = k-1

            dr = 0d0
            do j=1,manz
                dr = dr + rvec(j)*rvec(j)
            end do

            if (k.eq.1) then
                dr0  = dr*eps
                beta = 0d0
            else
                beta = dr/dr1
            end if

            if (dr.le.dr0) goto 10

            do j=1,manz
                pvec(j) = rvec(j) + beta*pvec(j)
            end do
            
            IF (ltri==0) THEN
               CALL bpdc(bvecdc,pvec)
            ELSE IF (ltri==1) THEN
               call bpdctri(bvecdc,pvec)
            END IF

            dr1 = 0d0
            do j=1,manz
                dr1 = dr1 + pvec(j)*bvecdc(j)
            end do

            alpha = dr/dr1

            do j=1,manz
                dpar(j) = dpar(j) + dcmplx(alpha*pvec(j))
                rvec(j) = rvec(j) - alpha*bvecdc(j)
            end do

            dr1 = dr

c Residuum speichern
            cgres(k+1) = real(eps*dr/dr0)
        end do

        ncg = ncgmax

c Anzahl an CG-steps speichern
10      cgres(1) = real(ncg)

        return
        end

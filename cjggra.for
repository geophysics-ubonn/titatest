        subroutine cjggra(bvec)

c Unterprogramm berechnet Modellverbesserung mittels konjugierter
c Gradienten.

c Andreas Kemna                                            01-Mar-1996
c                                       Letzte Aenderung   04-Feb-1998
         
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
        complex         * 16    rvec(mmax),pvec(mmax)

c Skalare
        complex         * 16    beta
        real            * 8     alpha,
     1                          dr,dr0,dr1

c Hilfsvariablen
        integer         * 4     j,k

c.....................................................................

        do j=1,manz
            dpar(j) = dcmplx(0d0)
            rvec(j) = bvec(j)
            pvec(j) = dcmplx(0d0)
        end do

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

c Fletcher-Reeves-Version
                beta = dcmplx(dr/dr1)

cakc Polak-Ribiere-Version
cak                beta = dcmplx(0d0)
cak                do j=1,manz
cak                    beta = beta + dconjg(bvec(j))*rvec(j)
cak                end do
cak                beta = beta*dcmplx(-alpha/dr1)
            end if

            do j=1,manz
                pvec(j) = rvec(j) + beta*pvec(j)
            end do

            call bp(bvec,pvec)

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

c Residuum speichern
            cgres(k+1) = real(eps*dr/dr0)
        end do

        ncg = ncgmax

c Anzahl an CG-steps speichern
10      cgres(1) = real(ncg)

        return
        end

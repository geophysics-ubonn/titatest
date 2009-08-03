        subroutine bpsto(bvec,pvec)

c Unterprogramm berechnet b = B * p . Angepasst an die neue Regularisierungsmatrix (stoch. Kovarianzmatrix).

c Andreas Kemna                                            29-Feb-1996
c                                       Letzte Aenderung   30-Jun-2009

c.....................................................................

        USE alloci

        INCLUDE 'parmax.fin'
        INCLUDE 'dat.fin'
        INCLUDE 'model.fin'
        INCLUDE 'fem.fin'
        INCLUDE 'inv.fin'
        INCLUDE 'konv.fin'

c.....................................................................

c EIN-/AUSGABEPARAMETER:

c Vektoren
        complex         * 16    bvec(mmax)
        complex         * 16    pvec(mmax)
        complex         * 16     pvecdum(manz)

c.....................................................................

c PROGRAMMINTERNE PARAMETER:

c Hilfsvektor
        complex         * 16    ap(nmax)

c Hilfsvariablen
        complex         * 16    cdum
        integer         * 4     i,j

c.....................................................................

c A * p  berechnen (skaliert)
        do i=1,nanz
            ap(i) = dcmplx(0d0)

            do j=1,manz
                ap(i) = ap(i) + pvec(j)*sens(i,j)*dcmplx(fak(j))
            end do
        end do

ccoaa R^m * p  berechnen (skaliert)

	do j=1,manz
           pvecdum(i)=pvec(i)*dcmplx(fak(i))
	end do

        bvec(1:manz)=MATMUL(dcmplx(smatm),pvecdum(1:manz))  

c A^h * R^d * A * p + l * R^m * p  berechnen (skaliert)
        do j=1,manz
            cdum = dcmplx(0d0)

            do i=1,nanz
                cdum = cdum + dconjg(sens(i,j))*
     1                        dcmplx(wmatd(i)*dble(wdfak(i)))*ap(i)
            end do

            bvec(j) = cdum + dcmplx(lam)*bvec(j)
            bvec(j) = bvec(j)*dcmplx(fak(j))
        end do

        return
        end

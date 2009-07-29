        subroutine bp(bvec,pvec)

c Unterprogramm berechnet b = B * p .

c Andreas Kemna                                            29-Feb-1996
c                                       Letzte Aenderung   10-Nov-1997

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

c R^m * p  berechnen (skaliert)
        do i=1,manz
            cdum = dcmplx(0d0)

            if (i.gt.1)
     1          cdum = pvec(i-1)*dcmplx(smatm(i-1,2)*fak(i-1))
            if (i.lt.manz)
     1          cdum = cdum + pvec(i+1)*dcmplx(smatm(i,2)*fak(i+1))
            if (i.gt.nx)
     1          cdum = cdum + pvec(i-nx)*dcmplx(smatm(i-nx,3)*fak(i-nx))
            if (i.lt.manz-nx+1)
     1          cdum = cdum + pvec(i+nx)*dcmplx(smatm(i,3)*fak(i+nx))

            bvec(i) = cdum + pvec(i)*dcmplx(smatm(i,1)*fak(i))
        end do

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

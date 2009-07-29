        subroutine bpdctri(bvec,pvec)

! Unterprogramm berechnet b = B * p .

! Andreas Kemna                                            29-Feb-1996
!                                       Letzte Aenderung   09-Jan-1998
! geaendert auf zunaechst reine Triangulation
! Roland Blaschek, 5.6.2003
! jetzt allgemein gueltig, 12.6.2003
!.....................................................................

        USE alloci

        INCLUDE 'parmax.fin'
        INCLUDE 'dat.fin'
        INCLUDE 'model.fin'
        INCLUDE 'fem.fin'
        INCLUDE 'inv.fin'
        INCLUDE 'konv.fin'

!.....................................................................

! EIN-/AUSGABEPARAMETER:

! Vektoren
        real            * 8     bvec(mmax)
        real            * 8     pvec(mmax)

!.....................................................................

! PROGRAMMINTERNE PARAMETER:

! Hilfsvektor
        real            * 8     ap(nmax)

! Hilfsvariablen
        real            * 8     dum
        integer         * 4     i,j

!.....................................................................
!WRITE(*,*) "BPDCtri",errnr 
! A * p  berechnen (skaliert)
        do i=1,nanz
            ap(i) = 0d0

            if (ldc) then
                do j=1,manz
                    ap(i) = ap(i) + pvec(j)*sensdc(i,j)*fak(j)
                end do
            else if (lip) then
                do j=1,manz
                    ap(i) = ap(i) + pvec(j)*dble(sens(i,j))*fak(j)
                end do
            end if
        end do

! R^m * p  berechnen (skaliert)
        DO i=1,manz
            dum = 0d0
            if ((.not.lfcon).OR.(.not.lme)) then
            DO j=1,nachbar(i,0)
            IF (nachbar(i,j)/=0) dum = dum + pvec(nachbar(i,j)) * 
     1                       smatm(i,j)*fak(nachbar(i,j))
            END DO
            else
            end if
            if (lfregu.OR.lfregud) then
            bvec(i) = dum + pvec(i)*(smatm(i,0)+
     1       1/(fvel**2 * fincr**2))*fak(i)
            
            elseif (lfcon) then
            bvec(i) = dum + pvec(i)*(
     1       1/(fvel**2 * fincr**2))*fak(i)
            elseif (lme) then
        bvec(i) = dum + pvec(i)*(1/exp(par(i)))*fak(i)
            elseif (lmesmooth) then
        bvec(i) = dum + pvec(i)*(smatm(i,0)+1/exp(par(i)))*fak(i)
            
            else
        bvec(i) = dum + pvec(i)*smatm(i,0)*fak(i)
        end if
        END DO

! A^h * R^d * A * p + l * R^m * p  berechnen (skaliert)
        do j=1,manz
            dum = 0d0

            if (ldc) then
                do i=1,nanz
                    dum = dum + sensdc(i,j)* 
     1                          wmatd(i)*dble(wdfak(i))*ap(i)
                end do
            else if (lip) then
                do i=1,nanz
                    dum = dum + dble(sens(i,j))* 
     1                           wmatd(i)*dble(wdfak(i))*ap(i)
                end do
            end if
            if (.not.lfregud) then
            bvec(j) = dum + lam*bvec(j)
            else 
            bvec(j) = cdum + dcmplx(lam)*bvec(j)+pvec(j)*(
     1   1/(fvel**2 * fincr**2))*fak(j)
            end if
            bvec(j) = bvec(j)*fak(j)
        end do

        return
        end

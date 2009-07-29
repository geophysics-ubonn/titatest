      subroutine bptri(bvec,pvec)

! Unterprogramm berechnet b = B * p .

! Andreas Kemna                                            29-Feb-1996
!                                       Letzte Aenderung   10-Nov-1997
! geaendert auf zunaechst reine Triangulation
! Roland Blaschek, 5.6.2003
! jetzt allgemein gueltig, 12.6.2003
!...................................................................

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
        complex         * 16    bvec(mmax)
        complex         * 16    pvec(mmax)

!.....................................................................

! PROGRAMMINTERNE PARAMETER:

! Hilfsvektor
        complex         * 16    ap(nmax)

! Hilfsvariablen
        complex         * 16    cdum
        integer         * 4     i,j

!.....................................................................
!
! A * p  berechnen (skaliert)
        do i=1,nanz
            ap(i) = dcmplx(0d0)

            do j=1,manz
                ap(i) = ap(i) + pvec(j)*sens(i,j)*dcmplx(fak(j))
            end do
        end do
        

! R^m * p  berechnen (skaliert)
      DO i=1,manz
       cdum = dcmplx(0d0)
       if ((.not.lfcon).OR.(.not.lme)) then
       DO j=1,nachbar(i,0)
       IF (nachbar(i,j)/=0) cdum = cdum + pvec(nachbar(i,j))* 
     1                          DCMPLX(smatm(i,j))*fak(nachbar(i,j))
       END DO
       else 
       end if
       if (lfregu.OR.lfregud) then

      bvec(i) = cdum + pvec(i)*(dcmplx(smatm(i,0))+ 
     1 1/(fvel**2 * fincr**2))*fak(i)
     
      elseif (lfcon) then
      bvec(i) = cdum + pvec(i)*(
     1 1/(fvel**2 * fincr**2))*fak(i)
        elseif (lme) then
        bvec(i) = cdum + pvec(i)*(1/exp(par(i)))*fak(i)
        elseif (lmesmooth) then
        bvec(i) = cdum + pvec(i)*(smatm(i,0)+1/exp(par(i)))*fak(i)
      else
      bvec(i) = cdum + pvec(i)*dcmplx(smatm(i,0))*fak(i)
        END IF
      END DO
    

! A^h * R^d * A * p + l * R^m * p  berechnen (skaliert)
        do j=1,manz
            cdum = dcmplx(0d0)

            do i=1,nanz
                cdum = cdum + dconjg(sens(i,j))* 
     1                         dcmplx(wmatd(i)*dble(wdfak(i)))*ap(i)
            end do
            if (.not.lfregud) then
            bvec(j) = cdum + dcmplx(lam)*bvec(j)
            else
            
            bvec(j) = cdum + dcmplx(lam)*bvec(j)+pvec(j)*(
     1   1/(fvel**2 * fincr**2))*fak(j)
            end if
            
            bvec(j) = bvec(j)*dcmplx(fak(j))
        end do
      

        return
        end


MODULE cg_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!$ This MODULE should deliver the interface for the Conjugate Gradient!  
!!!$ Method routines which are utilized to solve the normal equations   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!$ Copyright by Andreas Kemna 2010
!!!$
!!!$ Created by Roland Martin               30-Jul-2010
!!!$
!!!$ Last changed       RM                  Jul-2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE cjgmod
  USE alloci , ONLY : sens,sensdc,smatm,nachbar
  USE femmod , ONLY : fak,ldc
  USE elemmod, ONLY : smaxs
  USE invmod , ONLY : lip,wmatd,wdfak,dpar
  USE errmod , ONLY : errnr,fetxt
  USE konvmod , ONLY : ltri,lam,nx,nz
  USE modelmod , ONLY : manz
  USE datmod , ONLY : nanz

  IMPLICIT none

!!!$ DC subroutines
  PUBLIC :: cjggdc
!!!$ Subroutine calculates model update 
!!!$ with preconditioned conjugate gradient method
  PRIVATE :: bpdc
!!!$  subroutine calculates b = B * p (RHS) smooth regularization
  PRIVATE :: bpdctri
!!! same but for unstructured grids
  PRIVATE :: bpdclma
!!!$ for Levenberg and Levenberg-Marquardt damping
  PRIVATE :: bpdcsto
!!$ for stochastical regularization

!!$ IP subroutines
  PUBLIC :: cjggra
!!!$ Subroutine calculates model update for COMPLEX case
!!!$ with preconditioned conjugate gradient method
  PRIVATE :: bp
!!!$  subroutine calculates b = B * p (RHS) smooth regularization
  PRIVATE :: bptri
!!! same but for unstructured grids
  PRIVATE :: bplma
!!!$ for Levenberg and Levemnberg-Marquardt damping
  PRIVATE :: bpsto
!!$ for stochastical regularization, MATMUL is 
!!$ now explicitly formed because of conjugate complex


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                          DC_PART                               !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cjggdc()
!!!$    Unterprogramm berechnet Modellverbesserung mittels konjugierter
!!!$    Gradienten.
!!!$
!!!$    Andreas Kemna                                        01-Mar-1996
!!!$    Letzte Aenderung                                     29-Jul-2009
!!!$....................................................................
!!!$    PROGRAMMINTERNE PARAMETER:
!!!$
!!!$    Skalare
    REAL(KIND(0D0)) :: beta,alpha,dr,dr0,dr1
!!!$
!!!$    Hilfsvariablen
    INTEGER         :: j,k
!!!$
!!!$....................................................................
    ALLOCATE (rvecdc(manz),pvecdc(manz),apdc(nanz),&
         bvecdc(manz),stat=errnr)
    IF (errnr /= 0) THEN
       fetxt = 'Error memory allocation rve!!!$in cjggdc'
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

       WRITE (*,'(a,t40,I5,t55,G10.4,t70,G10.4)',ADVANCE='no')&
            ACHAR(13)//TRIM(fetxt),k,dr,dr0

       if (dr.le.dr0) goto 10

       pvecdc= rvecdc+ beta * pvecdc

       IF (ltri == 0) THEN
          CALL bpdc

       ELSE IF (ltri == 1.OR.ltri == 2.OR.&
            (ltri > 4 .AND. ltri < 15)) THEN
          CALL bpdctri

       ELSE IF (ltri == 3.OR.ltri == 4) THEN
          CALL bpdclma

       ELSE IF (ltri == 15) THEN
          CALL bpdcsto

       END IF

       dr1 = DOT_PRODUCT(pvecdc,bvecdc)

       alpha = dr/dr1

       dpar = dpar + DCMPLX(alpha) * DCMPLX(pvecdc)
       rvecdc= rvecdc- alpha * bvecdc

!!!$rm update speichern
       dr1 = dr

!!!$    Residuum speichern
       cgres(k+1) = real(eps*dr/dr0)
    end do

    ncg = ncgmax

!!!$    Anzahl an CG-steps speichern
10  cgres(1) = real(ncg)

    DEALLOCATE (pvecdc,rvecdc,apdc,bvecdc)

  end subroutine cjggdc

  subroutine bpdc()
!!$
!!!$    Unterprogramm berechnet b = B * p .
!!!$
!!!$    Andreas Kemna                                      29-Feb-1996
!!$
!!!$    Last changes   RM                                  Jul-2010
!!$
!!!$...................................................................
!!!$    PROGRAMMINTERNE PARAMETER:
!!!$    Hilfsvariablen
    REAL(KIND(0D0))    ::     dum
    INTEGER         ::     i,j

!!!$....................................................................

!!!$    A * p  berechnen (skaliert)
    do i=1,nanz
       apdc(i) = 0d0

       if (ldc) then
          do j=1,manz
             apdc(i) = apdc(i) + pvecdc(j)*sensdc(i,j)*fak(j)
          end do
       else if (lip) then
          do j=1,manz
             apdc(i) = apdc(i) + pvecdc(j)*dble(sens(i,j))*fak(j)
          end do
       end if
    end do

!!!$    R^m * p  berechnen (skaliert)
    do i=1,manz
       dum = 0d0

       if (i.gt.1) &
            dum = pvecdc(i-1)*smatm(i-1,2)*fak(i-1)
       if (i.lt.manz) &
            dum = dum + pvecdc(i+1)*smatm(i,2)*fak(i+1)
       if (i.gt.nx) &
            dum = dum + pvecdc(i-nx)*smatm(i-nx,3)*fak(i-nx)
       if (i.lt.manz-nx+1) &
            dum = dum + pvecdc(i+nx)*smatm(i,3)*fak(i+nx)

       bvecdc(i) = dum + pvecdc(i)*smatm(i,1)*fak(i)
    end do

!!!$    A^h * R^d * A * p + l * R^m * p  berechnen (skaliert)
    do j=1,manz
       dum = 0d0

       if (ldc) then
          do i=1,nanz
             dum = dum + sensdc(i,j) * &
                  wmatd(i)*dble(wdfak(i))*apdc(i)
          end do
       else if (lip) then
          do i=1,nanz
             dum = dum + dble(sens(i,j)) * &
                  wmatd(i)*dble(wdfak(i))*apdc(i)
          end do
       end if

       bvecdc(j) = dum + lam*bvecdc(j)
       bvecdc(j) = bvecdc(j)*fak(j)
    end do

  end subroutine bpdc

  subroutine bpdctri()
!!!$    
!!!$    Unterprogramm berechnet b = B * p .
!!!$    Fuer beliebige Triangulierung
!!!$    
!!!$    Copyright by Andreas Kemna         2009
!!!$
!!!$    Created by Roland Martin                            29-Jul-2009
!!!$     
!!!$    Last changes      RM                                   Jul-2010
!!!$    
!!!$....................................................................
!!!$.....................................................................
!!!$
!!!$!     PROGRAMMINTERNE PARAMETER:

!!!$!     Hilfsvariablen
    REAL(KIND(0D0))    ::     dum
    INTEGER         ::     i,j
!!!$!.....................................................................
!!!!     WRITE(*,*) "BPDCtri",errnr 
!!!!     A * p  berechnen (skaliert)

    do i=1,nanz
       apdc(i) = 0d0

       if (ldc) then
          do j=1,manz
             apdc(i) = apdc(i) + pvecdc(j) * sensdc(i,j) * fak(j)
          end do
       else if (lip) then
          do j=1,manz
             apdc(i) = apdc(i) + pvecdc(j) * dble(sens(i,j)) * fak(j)
          end do
       end if
    end do

    !     R^m * p  berechnen (skaliert)
    DO i=1,manz
       dum = 0d0
       DO j=1,smaxs
          IF (nachbar(i,j) /= 0) dum = dum + pvecdc(nachbar(i,j)) * & 
               smatm(i,j) * fak(nachbar(i,j)) ! off diagonals
       END DO
       !     main diagonal
       bvecdc(i) = dum + pvecdc(i) * smatm(i,smaxs+1) * fak(i) 
    END DO

    !     A^h * R^d * A * p + l * R^m * p  berechnen (skaliert)
    do j=1,manz
       dum = 0d0

       if (ldc) then
          do i=1,nanz
             dum = dum + sensdc(i,j) * wmatd(i) * dble(wdfak(i)) * &
                  apdc(i)
          end do
       else if (lip) then
          do i=1,nanz
             dum = dum + dble(sens(i,j)) * wmatd(i) * &
                  dble(wdfak(i)) * apdc(i)
          end do
       end if

       bvecdc(j) = dum + lam * bvecdc(j)

       bvecdc(j) = bvecdc(j) * fak(j)
    end do

  end subroutine bpdctri



  subroutine bpdclma()
!!!$    
!!!$    Unterprogramm berechnet b = B * p . 
!!!$    Angepasst an Levenberg-Marquardt-Daempfung
!!!$   
!!!$    Copyright by Andreas Kemna 2010
!!!$    
!!!$    Created by Roland Martin                            24-Feb-2010
!!!$    
!!!$    Last changes        RM                                Jul-2010
!!!$
!!!$....................................................................
!!!$    PROGRAMMINTERNE PARAMETER:

!!!$    Hilfsvariablen
    REAL(KIND(0D0))    ::     dum
    INTEGER         ::     i,j

!!!$....................................................................

!!!$    A * p  berechnen (skaliert)
    do i=1,nanz
       apdc(i) = 0d0

       if (ldc) then
          do j=1,manz
             apdc(i) = apdc(i) + pvecdc(j)*sensdc(i,j)*fak(j)
          end do
       else if (lip) then
          do j=1,manz
             apdc(i) = apdc(i) + pvecdc(j)*dble(sens(i,j))*fak(j)
          end do
       end if
    end do

!!!$    R^m * p  berechnen (skaliert)
    do i=1,manz
       bvecdc(i)=pvecdc(i)*fak(i)*smatm(i,1) ! damping stuff..
    end do

!!!$    A^h * R^d * A * p + l * R^m * p  berechnen (skaliert)
    do j=1,manz
       dum = 0d0

       if (ldc) then
          do i=1,nanz
             dum = dum + sensdc(i,j) * &
                  wmatd(i)*dble(wdfak(i))*apdc(i)
          end do
       else if (lip) then
          do i=1,nanz
             dum = dum + dble(sens(i,j)) * &
                  wmatd(i)*dble(wdfak(i))*apdc(i)
          end do
       end if

       bvecdc(j) = dum + lam*bvecdc(j)
       bvecdc(j) = bvecdc(j)*fak(j)
    end do

  end subroutine bpdclma

  subroutine bpdcsto()
!!!$
!!!$    Unterprogramm berechnet b = B * p . 
!!!$    Angepasst an die neue Regularisierungsmatrix
!!!$    (stoch. Kovarianzmatrix)
!!!$
!!!$    Copyright by Andreas Kemna 2009
!!!$    
!!!$    Created by Roland Martin                             10-Jun-2009
!!!$
!!!$    Last changes        RM                                Jul-2010
!!!$
!!!$....................................................................
!!!$    PROGRAMMINTERNE PARAMETER:
!!!$    Hilfsvariablen
    REAL(KIND(0D0)),ALLOCATABLE,DIMENSION(:) :: pvec2
    REAL(KIND(0D0))    ::     dum
    INTEGER         ::     i,j

!!!$....................................................................
!!!$      ALLOCATE (pvec2(manz),stat=errnr)
!!!$      IF (errnr /= 0) THEN
!!!$         fetxt = 'Error memory allocation pvec2 in bpdcsto'
!!!$         errnr = 94
!!!$         RETURN
!!!$      END IF

!!!$    A * p  berechnen (skaliert)
    do i=1,nanz
       apdc(i) = 0d0

       if (ldc) then
          do j=1,manz
             apdc(i) = apdc(i) + pvecdc(j)*sensdc(i,j)*fak(j)
          end do
       else if (lip) then
          do j=1,manz
             apdc(i) = apdc(i) + pvecdc(j)*dble(sens(i,j))*fak(j)
          end do
       end if
    end do

!!!$    R^m * p  berechnen (skaliert)
!!!$caa   Abgeändert auf (4 Zeilen)
!!!$      do i=1,manz
!!!$         pvec2(i)=pvecdc(i)*fak(i)
!!!$      end do
    do j = 1 , manz
       bvecdc(j) = 0.
       DO i = 1 , manz
          bvecdc(j) = bvecdc(j) + pvecdc(i) * smatm(i,j) * fak(i)
       END DO
    end do

!!!$
!!!$      bvecdc= MATMUL(smatm,pvec2)
!!!$    A^h * R^d * A * p + l * R^m * p  berechnen (skaliert)
    do j=1,manz
       dum = 0d0

       if (ldc) then
          do i=1,nanz
             dum = dum + sensdc(i,j) * &
                  wmatd(i)*dble(wdfak(i))*apdc(i)
          end do
       else if (lip) then
          do i=1,nanz
             dum = dum + dble(sens(i,j)) * &
                  wmatd(i)*dble(wdfak(i))*apdc(i)
          end do
       end if

       bvecdc(j) = dum + lam*bvecdc(j)
       bvecdc(j) = bvecdc(j)*fak(j)
    end do

    IF (ALLOCATED (pvec2)) DEALLOCATE (pvec2)

  end subroutine bpdcsto


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                         IP_PART                                !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cjggra()
!!!$    Unterprogramm berechnet Modellverbesserung mittels konjugierter
!!!$    Gradienten.
!!!$
!!!$    Andreas Kemna                                        01-Mar-1996
!!!$    Last changes   RM                                    Jul-2010
!!!$....................................................................
!!!$    PROGRAMMINTERNE PARAMETER:
!!!$    Skalare
    COMPLEX(KIND(0D0)) :: beta
    REAL(KIND(0D0))    :: alpha,dr,dr0,dr1
!!$
!!!$    Hilfsvariablen
    INTEGER            :: j,k
!!!$....................................................................

    ALLOCATE (rvec(manz),pvec(manz),ap(nanz),stat=errnr)
    IF (errnr /= 0) THEN
       fetxt = 'Error memory allocation rve!!!$in cjggdc'
       errnr = 94
       RETURN
    END IF

    dpar = 0.
    rvec = bvec
    pvec = 0.

    fetxt = 'CG iteration'

    do k=1,ncgmax

       ncg = k-1

!!!$         dr = 0d0
!!!$         do j=1,manz
!!!$            dr = dr + dble(dconjg(rvec(j))*rvec(j))
!!!$         end do

       dr = DOT_PRODUCT(DCONJG(rvec),rvec)

       if (k.eq.1) then
          dr0  = dr*eps
          beta = dcmplx(0d0)
       else
!!!$    Fletcher-Reeves-Version
          beta = dcmplx(dr/dr1)
!!!$    ak!!!$Polak-Ribiere-Version
!!!$    ak                beta = dcmplx(0d0)
!!!$    ak                do j=1,manz
!!!$    ak                    beta = beta + dconjg(bvec(j))*rvec(j)
!!!$    ak                end do
!!!$    ak                beta = beta*dcmplx(-alpha/dr1)
       END IF

       WRITE (*,'(a,t40,I5,t55,G10.4,t70,G10.4)',ADVANCE='no')&
            ACHAR(13)//TRIM(fetxt),k,dr,dr0

       if (dr.le.dr0) goto 10

       pvec = rvec + beta * pvec

       IF (ltri == 0) THEN
          CALL bp
       ELSE IF (ltri == 1.OR.ltri == 2.OR.&
            (ltri > 4 .AND. ltri < 15)) THEN
          call bptri
       ELSE IF (ltri == 3.OR.ltri == 4) THEN
          CALL bplma
       ELSE IF (ltri == 15) THEN
          call bpsto
       END IF

       dr1 = DOT_PRODUCT(DCONJG(pvec),bvec)
!!!$
!!!$         dr1 = 0d0
!!!$         do j=1,manz
!!!$            dr1 = dr1 + dble(dconjg(pvec(j))*bvec(j))
!!!$         end do

       alpha = dr/dr1

       dpar = dpar + DCMPLX(alpha) * pvec
       rvec = rvec - DCMPLX(alpha) * bvec

       dr1 = dr

!!!$    Residuum speichern
       cgres(k+1) = real(eps*dr/dr0)
    end do

    ncg = ncgmax

!!!$    Anzahl an CG-steps speichern
10  cgres(1) = real(ncg)

    DEALLOCATE (rvec,pvec,ap)

  end subroutine cjggra

  subroutine bp()
!!!$
!!!$    Unterprogramm berechnet b = B * p .
!!!$
!!!$    Andreas Kemna                                        29-Feb-1996
!!!$     
!!!$    Last changes      RM                                   Jul-2010
!!!$    
!!!$....................................................................
!!!$    PROGRAMMINTERNE PARAMETER:

!!!$    Hilfsvariablen
    COMPLEX(KIND(0D0)) ::    cdum
    INTEGER         ::     i,j

!!!$....................................................................

!!!$    A * p  berechnen (skaliert)
    do i=1,nanz
       ap(i) = dcmplx(0d0)

       do j=1,manz
          ap(i) = ap(i) + pvec(j)*sens(i,j)*dcmplx(fak(j))
       end do
    end do

!!!$    R^m * p  berechnen (skaliert)
    do i=1,manz
       cdum = dcmplx(0d0)

       if (i.gt.1) &
            cdum = pvec(i-1)*dcmplx(smatm(i-1,2)*fak(i-1))
       if (i.lt.manz) &
            cdum = cdum + pvec(i+1)*dcmplx(smatm(i,2)*fak(i+1))
       if (i.gt.nx) &
            cdum = cdum + pvec(i-nx)*dcmplx(smatm(i-nx,3)*fak(i-nx))
       if (i.lt.manz-nx+1) &
            cdum = cdum + pvec(i+nx)*dcmplx(smatm(i,3)*fak(i+nx))
       bvec(i) = cdum + pvec(i)*dcmplx(smatm(i,1)*fak(i))
    end do

!!!$    A^h * R^d * A * p + l * R^m * p  berechnen (skaliert)
    do j=1,manz
       cdum = dcmplx(0d0)

       do i=1,nanz
          cdum = cdum + dconjg(sens(i,j)) * &
               dcmplx(wmatd(i)*dble(wdfak(i)))*ap(i)
       end do

       bvec(j) = cdum + dcmplx(lam)*bvec(j)
       bvec(j) = bvec(j)*dcmplx(fak(j))
    end do

  end subroutine bp

  subroutine bptri()
!!!$    
!!!$    Unterprogramm berechnet b = B * p .
!!!$    Fuer beliebige Triangulierung
!!!$    
!!!$    Copyright by Andreas Kemna         2009
!!!$
!!!$    Created by Roland Martin                           29-Jul-2009
!!!$    
!!!$    Last changes      RM                                  Jul-2010
!!!$    
!!!$..................................................................
!!!$.....................................................................
!!!$
!!!$     PROGRAMMINTERNE PARAMETER:
!!!$
!!!$     Hilfsvariablen
    COMPLEX(KIND(0D0)) ::    cdum
    INTEGER         ::     i,j,idum
!!!.....................................................................
!!!     
!!!     A * p  berechnen (skaliert)

    do i=1,nanz
       ap(i) = dcmplx(0d0)

       do j=1,manz
          ap(i) = ap(i) + pvec(j)*sens(i,j)*dcmplx(fak(j))
       end do
    end do


    !     R^m * p  berechnen (skaliert)
    DO i=1,manz
       cdum = dcmplx(0d0)
       DO j=1,smaxs
          idum=nachbar(i,j)
          IF (idum/=0) cdum = cdum + pvec(idum) * & 
               DCMPLX(smatm(i,j)) * DCMPLX(fak(idum)) ! off diagonals
       END DO

       bvec(i) = cdum + pvec(i) * DCMPLX(smatm(i,smaxs+1)) * &
            DCMPLX(fak(i)) ! + main diagonal

    END DO


    !     A^h * R^d * A * p + l * R^m * p  berechnen (skaliert)
    do j=1,manz
       cdum = dcmplx(0d0)

       do i=1,nanz
          cdum = cdum + dconjg(sens(i,j))*dcmplx(wmatd(i) * &
               dble(wdfak(i)))*ap(i)
       end do

       bvec(j) = cdum + dcmplx(lam)*bvec(j)

       bvec(j) = bvec(j)*dcmplx(fak(j))
    end do

  end subroutine bptri

  subroutine bplma()
!!!$
!!!$    Unterprogramm berechnet b = B * p . 
!!!$    Angepasst an Levenberg-Marquardt-Daempfung
!!!$
!!!$    Copyright by Andreas Kemna       2010
!!!$    
!!!$    Created by Roland Martin                              24-Feb-2010
!!!$
!!!$    Last changes        RM                                Jul-2010
!!!$
!!!$....................................................................
!!!$    PROGRAMMINTERNE PARAMETER:
!!!$    Hilfsvariablen
    COMPLEX(KIND(0D0)) ::    cdum
    INTEGER         ::     i,j
!!!$....................................................................
!!!$    A * p  berechnen (skaliert)
    do i=1,nanz
       ap(i) = dcmplx(0d0)

       do j=1,manz
          ap(i) = ap(i) + pvec(j)*sens(i,j)*dcmplx(fak(j))
       end do
    end do

!!!$    coaa R^m * p  berechnen (skaliert)

    do j=1,manz
       bvec(i)=pvec(i)*dcmplx(fak(i))*DCMPLX(smatm(i,1))
    end do

!!!$    A^h * R^d * A * p + l * R^m * p  berechnen (skaliert)
    do j=1,manz
       cdum = dcmplx(0d0)

       do i=1,nanz
          cdum = cdum + dconjg(sens(i,j))* &
               dcmplx(wmatd(i)*dble(wdfak(i)))*ap(i)
       end do

       bvec(j) = cdum + dcmplx(lam)*bvec(j)
       bvec(j) = bvec(j)*dcmplx(fak(j))
    end do

  end subroutine bplma


  subroutine bpsto()
!!!$
!!!$    Unterprogramm berechnet b = B * p .
!!!$    Angepasst an die neue Regularisierungsmatrix 
!!!$    (stoch. Kovarianzmatrix) fuer komplexes Modell
!!!$
!!!$    Copyright by Andreas Kemna 2009
!!!$    
!!!$    Created by Roland Martin                              10-Jun-2009
!!!$
!!!$    Last changes   RM                                     Jul-2010
!!!$
!!!$....................................................................
!!!$    PROGRAMMINTERNE PARAMETER:
!!!$    Hilfsvariablen
    COMPLEX(KIND(0D0)) ::    cdum
    INTEGER         ::     i,j
!!!$....................................................................
!!!$    A * p  berechnen (skaliert)
    do i=1,nanz
       ap(i) = dcmplx(0d0)
       do j=1,manz
          ap(i) = ap(i) + pvec(j)*sens(i,j)*dcmplx(fak(j))
       end do
    end do

!!!$    R^m * p  berechnen (skaliert)
    DO j=1,manz
       bvec(j) = 0.
       DO i = 1,manz
          bvec(j) = bvec(j) + pvec(i) * DCMPLX(smatm(i,j)) * &
               DCMPLX(fak(j))
       END DO
    END DO

!!!$    A^h * R^d * A * p + l * R^m * p  berechnen (skaliert)
    do j=1,manz
       cdum = dcmplx(0d0)

       do i=1,nanz
          cdum = cdum + dconjg(sens(i,j)) * &
               dcmplx(wmatd(i)*dble(wdfak(i)))*ap(i)
       end do

       bvec(j) = cdum + dcmplx(lam)*bvec(j)
       bvec(j) = bvec(j)*dcmplx(fak(j))
    end do

  end subroutine bpsto


END MODULE cg_mod
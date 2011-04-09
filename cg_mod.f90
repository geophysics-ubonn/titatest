MODULE cg_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!$ This MODULE should deliver the interface for the Conjugate Gradient!  
!!!$ Method routines which are utilized to solve the normal equations   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!$ Copyright by Andreas Kemna 2010
!!!$
!!!$ Edited by Roland Martin               30-Jul-2010
!!!$
!!!$ Last changed       RM                  Feb-2011
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE alloci , ONLY : sens,sensdc,smatm,nachbar
  USE femmod , ONLY : fak,ldc
  USE elemmod, ONLY : smaxs
  USE invmod , ONLY : lip,wmatd,wdfak,dpar
  USE errmod , ONLY : errnr,fetxt
  USE konvmod , ONLY : ltri,lam,nx,nz,lverb
  USE modelmod , ONLY : manz
  USE datmod , ONLY : nanz
  USE cjgmod

  IMPLICIT none

  INTEGER,PARAMETER,PRIVATE :: ntd=2 ! number of threads
!!!$ we restrict ther thread numbers here to avoid atomization 
!!!$ of the problem since we are dealing with pure matrix vector product of types..

  PUBLIC :: cjg
!!!$ controls whather we have REAL or COMPLEX case
  

!!!$ DC subroutines
  PRIVATE :: cjggdc
!!!$ Subroutine calculates model update 
!!!$ with preconditioned conjugate gradient method

  PRIVATE :: bapdc
!!!$  sub calculates A * p (skaliert)
  PRIVATE :: bpdc
!!!$  subroutine calculates b = B * p (RHS) smooth regularization
  PRIVATE :: bpdctri
!!! same but for unstructured grids
  PRIVATE :: bpdclma
!!!$ for Levenberg and Levenberg-Marquardt damping
  PRIVATE :: bpdcsto
!!$ for stochastical regularization
  PRIVATE :: bbdc
!!!$ calculates  A^h * R^d * A * p + l * R^m * p  (skaliert)


!!$ IP subroutines
  PRIVATE :: cjggra
!!!$ Subroutine calculates model update for COMPLEX case
!!!$ with preconditioned conjugate gradient method

  PRIVATE :: bap
!!!$  sub calculates A * p (skaliert)
  PRIVATE :: bp
!!!$  subroutine calculates b = B * p (RHS) smooth regularization
  PRIVATE :: bptri
!!! same but for unstructured grids
  PRIVATE :: bplma
!!!$ for Levenberg and Levemnberg-Marquardt damping
  PRIVATE :: bpsto
!!$ for stochastical regularization, MATMUL is 
!!$ now explicitly formed because of conjugate complex
  PRIVATE :: bb
!!!$ calculates  A^h * R^d * A * p + l * R^m * p  (skaliert)


CONTAINS

  SUBROUTINE cjg
    if (ldc.or.lip) then
       CALL con_cjgmod (2,fetxt,errnr)
       IF (errnr /= 0) RETURN
       call cjggdc
       CALL des_cjgmod (2,fetxt,errnr)
       IF (errnr /= 0) RETURN
    else
       CALL con_cjgmod (3,fetxt,errnr)
       IF (errnr /= 0) RETURN
       call cjggra
       CALL des_cjgmod (3,fetxt,errnr)
       IF (errnr /= 0) RETURN
    end if

  END SUBROUTINE cjg
  

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
    INTEGER         :: k
!!!$
!!!$....................................................................

    if (lip) then
       bvecdc = dimag(bvec)
    else
       bvecdc = dble(bvec)
    end if

    dpar = DCMPLX(0D0)
    rvecdc = bvecdc
    pvecdc = 0D0

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

       IF (lverb) WRITE (*,'(a,t40,I5,t55,G10.4,t70,G10.4)',&
            ADVANCE='no')ACHAR(13)//TRIM(fetxt),k,dr,dr0

       if (dr.le.dr0) goto 10

       pvecdc= rvecdc + beta * pvecdc

       CALL bapdc

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

       CALL bbdc

       dr1 = DOT_PRODUCT(pvecdc,bvecdc) ! this is ok for ERT

       alpha = dr/dr1

       dpar = dpar + DCMPLX(alpha) * DCMPLX(pvecdc)
       rvecdc= rvecdc - alpha * bvecdc

!!!$rm update speichern
       dr1 = dr

!!!$    Residuum speichern
       cgres(k+1) = real(eps*dr/dr0)

    end do

    ncg = ncgmax

!!!$    Anzahl an CG-steps speichern
10  cgres(1) = real(ncg)

!    DEALLOCATE (pvecdc,rvecdc,apdc,bvecdc)
    
    RETURN

  end subroutine cjggdc

  SUBROUTINE bapdc
!!$
!!!$    Unterprogramm berechnet Hilfsvektor A * p (skaliert).
!!!$
!!!$    Andreas Kemna                                      29-Feb-1996
!!$
!!!$    Last changes   RM                                  Mar-2011
!!$
!!!$...................................................................
!!!$    PROGRAMMINTERNE PARAMETER:
!!!$    Hilfsvariablen
    INTEGER         ::     i,j

!!!$....................................................................

    apdc = 0D0

    !!$!$OMP PARALLEL NUM_THREADS (ntd) DEFAULT(none) &
    !!$!$OMP SHARED (nanz,apdc,ldc,lip,manz,pvecdc,sensdc,cgfac,sens)
    !!$!$OMP DO
!!!$    A * p  berechnen (skaliert)
    do i=1,nanz
       if (ldc) then
          do j=1,manz
             apdc(i) = apdc(i) + pvecdc(j)*sensdc(i,j)*cgfac(j)
          end do
       else if (lip) then
          do j=1,manz
             apdc(i) = apdc(i) + pvecdc(j)*dble(sens(i,j))*cgfac(j)
          end do
       end if
    end do
    !!$!$OMP END PARALLEL
  END SUBROUTINE bapdc

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

!!!$    R^m * p  berechnen (skaliert)
    do i=1,manz
       dum = 0d0

       if (i.gt.1) &
            dum = pvecdc(i-1)*smatm(i-1,2)*cgfac(i-1)
       if (i.lt.manz) &
            dum = dum + pvecdc(i+1)*smatm(i,2)*cgfac(i+1)
       if (i.gt.nx) &
            dum = dum + pvecdc(i-nx)*smatm(i-nx,3)*cgfac(i-nx)
       if (i.lt.manz-nx+1) &
            dum = dum + pvecdc(i+nx)*smatm(i,3)*cgfac(i+nx)

       bvecdc(i) = dum + pvecdc(i)*smatm(i,1)*cgfac(i)
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

    !     R^m * p  berechnen (skaliert)
    DO i=1,manz
       dum = 0d0
       DO j=1,smaxs
          IF (nachbar(i,j) /= 0) dum = dum + pvecdc(nachbar(i,j)) * & 
               smatm(i,j) * cgfac(nachbar(i,j)) ! off diagonals
       END DO
       !     main diagonal
       bvecdc(i) = dum + pvecdc(i) * smatm(i,smaxs+1) * cgfac(i) 
    END DO

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


!!!$    R^m * p  berechnen (skaliert)
    do i=1,manz
       bvecdc(i)=pvecdc(i)*cgfac(i)*smatm(i,1) ! damping stuff..
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
    REAL(KIND(0D0))    ::     dum
    INTEGER         ::     i,j

    bvecdc = 0D0
!!!$    R^m * p  berechnen (skaliert)
    !!$!$OMP PARALLEL NUM_THREADS (ntd) DEFAULT(none) PRIVATE (i,dum) &
    !!$!$OMP SHARED (manz,bvecdc,pvecdc,cgfac,smatm)
    !!$!$OMP DO
    do j = 1 , manz
       DO i = j , manz
          dum = pvecdc(i) * smatm(i,j) * cgfac(i)
          IF (i == j) THEN
             bvecdc(j) = bvecdc(j) + dum
          ELSE
             bvecdc(j) = bvecdc(j) + 2D0 * dum
          END IF
       END DO
    end do
    !!$!$OMP END PARALLEL
  end subroutine bpdcsto

  SUBROUTINE bbdc
!!$
!!!$    Unterprogramm berechnet A^h * R^d * A * p + l * R^m * p  berechnen (skaliert)
!!!$
!!!$    Andreas Kemna                                      29-Feb-1996
!!$
!!!$    Last changes   RM                                  Mar-2011
!!$
!!!$...................................................................
!!!$    PROGRAMMINTERNE PARAMETER:
!!!$    Hilfsvariablen
    REAL(KIND(0D0))    ::     dum
    INTEGER         ::     i,j

!!!$....................................................................

    !!$!$OMP PARALLEL NUM_THREADS (ntd) DEFAULT(none) PRIVATE (dum) &
    !!$!$OMP SHARED (manz,ldc,lip,nanz,sensdc,wmatd,wdfak,apdc,sens,bvecdc,lam,cgfac)
    !!$!$OMP DO
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
       bvecdc(j) = bvecdc(j)*cgfac(j)
    end do
    !!$!$OMP END PARALLEL

  end subroutine bbdc


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
!    COMPLEX(KIND(0D0)) :: beta
    REAL(KIND(0D0))    :: alpha,dr,dr0,dr1,beta
!!$
!!!$    Hilfsvariablen
    INTEGER            :: k,j
!!!$....................................................................


    dpar = dcmplx(0d0)
    rvec = bvec
    pvec = dcmplx(0d0)

    fetxt = 'CG iteration'

    do k=1,ncgmax

       ncg = k-1

       dr = 0d0
       DO j=1,manz
          dr = dr + DBLE(DCONJG(rvec(j)) * rvec(j))
       END DO

       if (k.eq.1) then
          dr0  = dr*eps
          beta = 0d0
       else
          if (dr.le.dr0) goto 10
!!!$    Fletcher-Reeves-Version
          beta = dr/dr1
!!!$    ak!!!$Polak-Ribiere-Version
!!$          beta = 0d0
!!$          do j=1,manz
!!$             beta = beta + dconjg(bvec(j))*rvec(j)
!!$          end do
!!$          beta = beta * -alpha/dr1
       END IF

       IF (lverb) WRITE (*,'(a,t40,I5,t55,G10.4,t70,G10.4)',&
            ADVANCE='no')ACHAR(13)//TRIM(fetxt),k,dr,dr0

       pvec = rvec + DCMPLX(beta) * pvec

       CALL bap

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

       CALL bb

       dr1 = 0d0
       DO j=1,manz
          dr1 = dr1 + DBLE(DCONJG(pvec(j)) * bvec(j))
       END DO

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

    
  end subroutine cjggra

  SUBROUTINE bap
!!!$
!!!$    Unterprogramm berechnet A * p  berechnen (skaliert)
!!!$
!!!$    Andreas Kemna                                        29-Feb-1996
!!!$     
!!!$    Last changes      RM                                   Mar-2011
!!!$    
!!!$....................................................................
!!!$    PROGRAMMINTERNE PARAMETER:

!!!$    Hilfsvariablen
    COMPLEX(KIND(0D0)) ::    cdum
    INTEGER         ::     i,j

!!!$....................................................................

!!!$    A * p  berechnen (skaliert)
    !!$!$OMP PARALLEL NUM_THREADS (ntd) DEFAULT(none) &
    !!$!$OMP SHARED (nanz,ap,manz,pvec,sens,cgfac)
    !!$!$OMP DO
    do i=1,nanz
       ap(i) = dcmplx(0d0)
       
       do j=1,manz
          ap(i) = ap(i) + pvec(j)*sens(i,j)*dcmplx(cgfac(j))
       end do
    end do
    !!$!$OMP END PARALLEL

  END SUBROUTINE bap

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
!!!$    R^m * p  berechnen (skaliert)
    do i=1,manz
       cdum = dcmplx(0d0)

       if (i.gt.1) &
            cdum = pvec(i-1)*dcmplx(smatm(i-1,2)*cgfac(i-1))
       if (i.lt.manz) &
            cdum = cdum + pvec(i+1)*dcmplx(smatm(i,2)*cgfac(i+1))
       if (i.gt.nx) &
            cdum = cdum + pvec(i-nx)*dcmplx(smatm(i-nx,3)*cgfac(i-nx))
       if (i.lt.manz-nx+1) &
            cdum = cdum + pvec(i+nx)*dcmplx(smatm(i,3)*cgfac(i+nx))
       bvec(i) = cdum + pvec(i)*dcmplx(smatm(i,1)*cgfac(i))
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
    !     R^m * p  berechnen (skaliert)
    DO i=1,manz
       cdum = dcmplx(0d0)
       DO j=1,smaxs
          idum=nachbar(i,j)
          IF (idum/=0) cdum = cdum + pvec(idum) * & 
               DCMPLX(smatm(i,j)) * cgfac(idum) ! off diagonals
       END DO

       bvec(i) = cdum + pvec(i) * DCMPLX(smatm(i,smaxs+1)) * &
            cgfac(i) ! + main diagonal

    END DO

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
!!!$    coaa R^m * p  berechnen (skaliert)

    bvec = pvec * DCMPLX(cgfac * smatm(:,1))
!!$    do j=1,manz
!!$       bvec(i)=pvec(i)*dcmplx(cgfac(i))*DCMPLX(smatm(i,1))
!!$    end do

  end subroutine bplma


  subroutine bpsto()
!!!$
!!!$    Unterprogramm berechnet b = B * p .
!!!$    Angepasst an die neue Regularisierungsmatrix 
!!!$    (stoch. Kovarianzmatrix) fuer komplexes Modell
!!!$
!!!$   TODO:
!!!$      since smatm is symmetric, it would be good to 
!!!$      exploit this..
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
!!!$    R^m * p  berechnen (skaliert)
    !!$!$OMP PARALLEL NUM_THREADS (ntd) DEFAULT(none) PRIVATE (i,cdum) &
    !!$!$OMP SHARED (manz,bvec,pvec,cgfac,smatm)
    !!$!$OMP DO
    DO j=1, manz
       bvec(j) = DCMPLX(0D0)
       DO i = j, manz
          cdum = pvec(i) * DCMPLX(smatm(i,j)) * DCMPLX(cgfac(j))
          IF (i == j) THEN
             bvec(j) = bvec(j) + cdum
          ELSE
             bvec(j) = bvec(j) + 2D0 * cdum
          END IF
       END DO
    END DO
    !!$!$OMP END PARALLEL

  end subroutine bpsto

  SUBROUTINE bb
!!!$
!!!$    Unterprogramm berechnet A^h * R^d * A * p + l * R^m * p  berechnen (skaliert)
!!!$
!!!$    Andreas Kemna                                        29-Feb-1996
!!!$     
!!!$    Last changes      RM                                   Mar-2011
!!!$    
!!!$....................................................................
!!!$    PROGRAMMINTERNE PARAMETER:

!!!$    Hilfsvariablen
    COMPLEX(KIND(0D0)) ::    cdum
    INTEGER         ::     i,j

!!!$....................................................................
!!!$    
    !!$!$OMP PARALLEL NUM_THREADS (ntd) DEFAULT(none) PRIVATE (cdum) &
    !!$!$OMP SHARED (manz,nanz,sens,wmatd,wdfak,ap,bvec,lam,cgfac)
    !!$!$OMP DO

    do j=1,manz
       cdum = dcmplx(0d0)

       do i=1,nanz
          cdum = cdum + dconjg(sens(i,j)) * &
               dcmplx(wmatd(i)*dble(wdfak(i)))*ap(i)
       end do

       bvec(j) = cdum + dcmplx(lam)*bvec(j)
       bvec(j) = bvec(j)*dcmplx(cgfac(j))
    end do

    !!$!$OMP END PARALLEL

  END SUBROUTINE bb


END MODULE cg_mod

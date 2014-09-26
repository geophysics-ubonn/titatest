MODULE cg_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This MODULE should deliver the interface for the Conjugate Gradient  
!! Method routines which are utilized to solve the normal equations   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!$ Copyright by Andreas Kemna 2010
!!!$
!!!$ Edited by Roland Martin               30-Jul-2010
!!!$
!!!$ Last changed       RM                  Feb-2011
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE alloci , ONLY : sens,sensdc,smatm,prec,lverbose
  USE femmod , ONLY : fak,ldc
  USE elemmod, ONLY : smaxs,nachbar
  USE invmod , ONLY : lfpi,wmatd,wdfak,dpar
  USE errmod , ONLY : errnr,fetxt
  USE konvmod , ONLY : ltri,lam,nx,nz,lverb,lw_ref, lam_ref, lam_ref_sw
  USE modelmod , ONLY : manz,w_ref_re,w_ref_im
  USE datmod , ONLY : nanz
  USE cjgmod

  IMPLICIT NONE

!> number of threads.
!! NOTE that we restrict the total thread numbers here to avoid atomization 
!! of the problem since we are dealing with pure matrix vector product
  INTEGER,PARAMETER,PRIVATE :: ntd=2 

!> controls whather we have REAL or COMPLEX case
  PUBLIC :: cjg
!> Subroutine calculates model update (DC)
!! with preconditioned conjugate gradient method
  PRIVATE :: cjggdc
!!  sub calculates A * p (skaliert)
  PRIVATE :: bapdc
!! calculates  A^h * R^d * A * p + l * R^m * p  (skaliert)
  PRIVATE :: bbdc
!! calculates b = B * p (RHS) smooth regularization
  PRIVATE :: bpdc
!! calculates b = B * p (RHS) same but for unstructured grids
  PRIVATE :: bpdctri
!! calculates b = B * p (RHS) for Levenberg and Levenberg-Marquardt damping
  PRIVATE :: bpdclma
!!$ calculates b = B * p (RHS) for stochastical regularization
  PRIVATE :: bpdcsto
!!$ adds additional parts to b = B * p (RHS) for reference model regularization
  PRIVATE :: bpdcref


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
!!$ adds additional parts to b = B * p (RHS) for reference model regularization
  PRIVATE :: bpref
  PRIVATE :: bb
!!!$ calculates  A^h * R^d * A * p + l * R^m * p  (skaliert)


CONTAINS
!> cjg flow control subroutine is called from outside
!! and checks for the different cases (DC,IP,FPI)
  SUBROUTINE cjg
    IF (ldc) THEN
       CALL con_cjgmod (2,fetxt,errnr)
       IF (errnr /= 0) RETURN
       CALL cjggdc
       CALL des_cjgmod (2,fetxt,errnr)
       IF (errnr /= 0) RETURN
    ELSE
    if (lverbose) print*,'calling cjggra'
       CALL con_cjgmod (3,fetxt,errnr)
       IF (errnr /= 0) RETURN
       CALL cjggra
       CALL des_cjgmod (3,fetxt,errnr)
       IF (errnr /= 0) RETURN
    END IF

  END SUBROUTINE cjg


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                          DC_PART                               !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE cjggdc()
!!!$    Unterprogramm berechnet Modellverbesserung mittels konjugierter
!!!$    Gradienten.
!!!$
!!!$    Andreas Kemna                                        01-Mar-1996
!!!$    Letzte Aenderung                                     29-Jul-2009
!!!$....................................................................
!!!$    PROGRAMMINTERNE PARAMETER:
!!!$
!!!$    Skalare
    REAL(prec) :: beta,alpha,dr,dr0,dr1
!!!$
!!!$    Hilfsvariablen
    INTEGER         :: k
!!!$
!!!$....................................................................

    IF (lfpi) THEN
       bvecdc = aimag(bvec)
    ELSE
       bvecdc = REAL(bvec)
    END IF

    dpar = dcmplx(0D0)
    rvecdc = bvecdc
    pvecdc = 0D0

    fetxt = 'CG iteration'
    DO k=1,ncgmax

       ncg = k-1

       dr = DOT_PRODUCT(rvecdc,rvecdc)

       IF (k.EQ.1) THEN
          dr0  = dr*eps
          beta = 0d0
       ELSE
          beta = dr/dr1
       END IF

       IF (lverb) WRITE (*,'(a,t40,I5,t55,G10.4,t70,G10.4)',&
            ADVANCE='no')ACHAR(13)//TRIM(fetxt),k,dr,dr0
       IF (dr.LE.dr0) GOTO 10

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
       
       IF (lw_ref) CALL bpdcref

       CALL bbdc

       dr1 = DOT_PRODUCT(pvecdc,bvecdc) ! this is ok for ERT

       alpha = dr/dr1

       dpar = dpar + dcmplx(alpha) * dcmplx(pvecdc)
       rvecdc= rvecdc - alpha * bvecdc

!!!$rm update speichern
       dr1 = dr

!!!$    Residuum speichern
       cgres(k+1) = REAL(eps*dr/dr0)

    END DO

    ncg = ncgmax

!!!$    Anzahl an CG-steps speichern
10  cgres(1) = REAL(ncg)

    !    DEALLOCATE (pvecdc,rvecdc,apdc,bvecdc)

    RETURN

  END SUBROUTINE cjggdc

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

    
    !$OMP PARALLEL NUM_THREADS (ntd) DEFAULT(none) &
    !$OMP SHARED (nanz,apdc,ldc,lfpi,manz,pvecdc,sensdc,cgfac,sens)

    apdc = 0D0

    !$OMP DO
!!!$    A * p  berechnen (skaliert)
    DO i=1,nanz
       IF (ldc) THEN
          DO j=1,manz
             apdc(i) = apdc(i) + pvecdc(j)*sensdc(i,j)*cgfac(j)
          END DO
       ELSE IF (lfpi) THEN
          DO j=1,manz
             apdc(i) = apdc(i) + pvecdc(j)*REAL(sens(i,j))*cgfac(j)
          END DO
       END IF
    END DO
    !$OMP END PARALLEL
  END SUBROUTINE bapdc

  SUBROUTINE bpdc()
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
    REAL(prec)    ::     dum
    INTEGER         ::     i,j

!!!$....................................................................

!!!$    R^m * p  berechnen (skaliert)
    DO i=1,manz
       dum = 0d0

       IF (i.GT.1) &
            dum = pvecdc(i-1)*smatm(i-1,2)*cgfac(i-1)
       IF (i.LT.manz) &
            dum = dum + pvecdc(i+1)*smatm(i,2)*cgfac(i+1)
       IF (i.GT.nx) &
            dum = dum + pvecdc(i-nx)*smatm(i-nx,3)*cgfac(i-nx)
       IF (i.LT.manz-nx+1) &
            dum = dum + pvecdc(i+nx)*smatm(i,3)*cgfac(i+nx)

       bvecdc(i) = dum + pvecdc(i)*smatm(i,1)*cgfac(i)
    END DO

  END SUBROUTINE bpdc

  SUBROUTINE bpdctri()
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
    REAL(prec)    ::     dum
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

  END SUBROUTINE bpdctri



  SUBROUTINE bpdclma()
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
    REAL(prec)    ::     dum
    INTEGER         ::     i,j

!!!$....................................................................


!!!$    R^m * p  berechnen (skaliert)
    DO i=1,manz
       bvecdc(i)=pvecdc(i)*cgfac(i)*smatm(i,1) ! damping stuff..
    END DO

  END SUBROUTINE bpdclma

  SUBROUTINE bpdcsto()
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
    REAL(prec)    ::     dum
    INTEGER         ::     i,j

!!!$    R^m * p  berechnen (skaliert)

!    !$OMP WORKSHARE
    bvecdc = MATMUL(smatm,pvecdc)
!    !$OMP END WORKSHARE

    bvecdc = bvecdc * cgfac

  END SUBROUTINE bpdcsto

  SUBROUTINE bpdcref()
!!!$    
!!!$    additional entries due to b = B*p
!!!$    for reference model regu
!!!$   
!!!$    Copyright by Andreas Kemna 2010
!!!$    
!!!$    Created by Roland Martin                            05-Sep-2012
!!!$    
!!!$    Last changes        RM                                Sep-2012
!!!$
!!!$....................................................................
!!!$    PROGRAMMINTERNE PARAMETER:

!!!$    Hilfsvariablen
    INTEGER         ::     i
    REAL (prec) :: rdum
!!!$....................................................................


!!!$    R^m * p  berechnen (skaliert)
    DO i=1,manz

       IF (w_ref_re(i) <= EPSILON(rdum)) CYCLE 

       rdum = pvecdc(i)*w_ref_re(i)

       bvecdc(i) = bvecdc(i) + rdum * lam_ref * cgfac(i)
!!!!$! according to damping stuff..
    END DO

  END SUBROUTINE bpdcref

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
    REAL(prec)    ::     dum
    INTEGER         ::     i,j

!!!$....................................................................

    !$OMP PARALLEL NUM_THREADS (ntd) DEFAULT(none) PRIVATE (dum) &
    !$OMP SHARED (manz,ldc,lfpi,nanz,sensdc,wmatd,wdfak,apdc,sens,bvecdc,lam,cgfac)
    !$OMP DO
    DO j=1,manz
       dum = 0d0

       IF (ldc) THEN
          DO i=1,nanz
             dum = dum + sensdc(i,j) * &
                  wmatd(i)*REAL(wdfak(i))*apdc(i)
          END DO
       ELSE IF (lfpi) THEN
          DO i=1,nanz
             dum = dum + REAL(sens(i,j)) * &
                  wmatd(i)*REAL(wdfak(i))*apdc(i)
          END DO
       END IF

       bvecdc(j) = dum + lam*bvecdc(j)
       bvecdc(j) = bvecdc(j)*cgfac(j)
    END DO
    !$OMP END PARALLEL

  END SUBROUTINE bbdc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                         IP_PART                                !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE cjggra()
!!!$    Unterprogramm berechnet Modellverbesserung mittels konjugierter
!!!$    Gradienten.
!!!$
!!!$    Andreas Kemna                                        01-Mar-1996
!!!$    Last changes   RM                                    Jul-2010
!!!$....................................................................
!!!$    PROGRAMMINTERNE PARAMETER:
!!!$    Skalare
!!$!!    COMPLEX(prec) :: beta
    REAL(prec)    :: alpha,dr,dr0,dr1,beta
!!$
!!!$    Hilfsvariablen
    INTEGER            :: k,j
!!!$....................................................................


    dpar = dcmplx(0d0)
    rvec = bvec
    pvec = dcmplx(0d0)

    fetxt = 'CG iteration'

    DO k=1,ncgmax

       ncg = k-1

       dr = 0d0
       DO j=1,manz
          dr = dr + REAL(dconjg(rvec(j)) * rvec(j))
       END DO

       IF (k.EQ.1) THEN
          dr0  = dr*eps
          beta = 0d0
       ELSE
        
!!!$    Fletcher-Reeves-Version
          beta = dr/dr1
!!!$    ak!!!$Polak-Ribiere-Version
!!$          beta = 0d0
!!$          do j=1,manz
!!$             beta = beta + dconjg(bvec(j))*rvec(j)
!!$          end do
!!$          beta = beta * -alpha/dr1
       END IF
!  IF (dr.LE.dr0) GOTO 10
       IF (lverb) WRITE (*,'(a,t40,I5,t55,G10.4,t70,G10.4)',&
            ADVANCE='no')ACHAR(13)//TRIM(fetxt),k,dr,dr0

       pvec = rvec + dcmplx(beta) * pvec
       CALL bap

       IF (ltri == 0) THEN
          CALL bp
       ELSE IF (ltri == 1.OR.ltri == 2.OR.&
            (ltri > 4 .AND. ltri < 15)) THEN
          CALL bptri
       ELSE IF (ltri == 3.OR.ltri == 4) THEN
          CALL bplma
       ELSE IF (ltri == 15) THEN
          CALL bpsto
       END IF
       IF (lw_ref) CALL bpref

       CALL bb
       dr1 = 0d0
       DO j=1,manz
          dr1 = dr1 + REAL(dconjg(pvec(j)) * bvec(j))
       END DO

       alpha = dr/dr1
       dpar = dpar + dcmplx(alpha) * pvec
       rvec = rvec - dcmplx(alpha) * bvec

       dr1 = dr

!!!$    Residuum speichern
       cgres(k+1) = REAL(eps*dr/dr0)

    END DO

    ncg = ncgmax

!!!$    Anzahl an CG-steps speichern
10  cgres(1) = REAL(ncg)


  END SUBROUTINE cjggra

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
    COMPLEX(prec) ::    cdum
    INTEGER         ::     i,j

!!!$....................................................................

!!!$    A * p  berechnen (skaliert)
    !$OMP PARALLEL NUM_THREADS (ntd) DEFAULT(none) &
    !$OMP SHARED (nanz,ap,manz,pvec,sens,cgfac)
    !$OMP DO
    DO i=1,nanz
       ap(i) = dcmplx(0d0)

       DO j=1,manz
          ap(i) = ap(i) + pvec(j)*sens(i,j)*dcmplx(cgfac(j))
       END DO
    END DO
    !$OMP END PARALLEL

  END SUBROUTINE bap

  SUBROUTINE bp()
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
    COMPLEX(prec) ::    cdum
    INTEGER         ::     i,j

!!!$....................................................................
!!!$    R^m * p  berechnen (skaliert)
    DO i=1,manz
       cdum = dcmplx(0d0)

       IF (i.GT.1) &
            cdum = pvec(i-1)*dcmplx(smatm(i-1,2)*cgfac(i-1))
       IF (i.LT.manz) &
            cdum = cdum + pvec(i+1)*dcmplx(smatm(i,2)*cgfac(i+1))
       IF (i.GT.nx) &
            cdum = cdum + pvec(i-nx)*dcmplx(smatm(i-nx,3)*cgfac(i-nx))
       IF (i.LT.manz-nx+1) &
            cdum = cdum + pvec(i+nx)*dcmplx(smatm(i,3)*cgfac(i+nx))
       bvec(i) = cdum + pvec(i)*dcmplx(smatm(i,1)*cgfac(i))
    END DO
  END SUBROUTINE bp

  SUBROUTINE bptri()
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
    COMPLEX(prec) ::    cdum
    INTEGER         ::     i,j,idum
!!!.....................................................................
    !     R^m * p  berechnen (skaliert)
    DO i=1,manz
       cdum = dcmplx(0d0)
       DO j=1,smaxs
          idum=nachbar(i,j)
          IF (idum/=0) cdum = cdum + pvec(idum) * & 
               dcmplx(smatm(i,j)) * cgfac(idum) ! off diagonals
       END DO

       bvec(i) = cdum + pvec(i) * dcmplx(smatm(i,smaxs+1)) * &
            cgfac(i) ! + main diagonal

    END DO

  END SUBROUTINE bptri

  SUBROUTINE bplma()
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
    COMPLEX(prec) ::    cdum
    INTEGER         ::     i,j
!!!$....................................................................
!!!$    coaa R^m * p  berechnen (skaliert)

    bvec = pvec * dcmplx(cgfac * smatm(:,1))
!!$    do j=1,manz
!!$       bvec(i)=pvec(i)*dcmplx(cgfac(i))*dcmplx(smatm(i,1))
!!$    end do

  END SUBROUTINE bplma


  SUBROUTINE bpsto()
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
    COMPLEX(prec) ::    cdum
    INTEGER         ::     i,j
!!!$....................................................................
!!!$    R^m * p  berechnen (skaliert)

    !$OMP WORKSHARE
    bvec = MATMUL(dcmplx(smatm),pvec)
    !$OMP END WORKSHARE

    bvec = bvec * dcmplx(cgfac)
  END SUBROUTINE bpsto

  SUBROUTINE bpref()
!!!$    
!!!$    additional entries due to b = B*p
!!!$    for reference model regu
!!!$   
!!!$    Copyright by Andreas Kemna 2010
!!!$    
!!!$    Created by Roland Martin                            05-Sep-2012
!!!$    
!!!$    Last changes        RM                                Sep-2012
!!!$
!!!$....................................................................
!!!$    PROGRAMMINTERNE PARAMETER:

!!!$    Hilfsvariablen
    INTEGER         ::     i
    COMPLEX (prec) :: cdum
!!!$....................................................................


!!!$    R^m * p  berechnen (skaliert)
    DO i=1,manz

       IF (w_ref_re(i) <= EPSILON(w_ref_re(i)) .AND. &
            w_ref_im(i) <= EPSILON(w_ref_im(i))) CYCLE 
!!$ scaling for real and aimaginary part separately       
       cdum = dcmplx(REAL(pvec(i))*w_ref_re(i), aimag(pvec(i))*w_ref_re(i))

       bvec(i) = bvec(i) + cdum * dcmplx(lam_ref * cgfac(i))

!!!!$! according to damping stuff..
    END DO

  END SUBROUTINE bpref

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
    COMPLEX(prec) ::    cdum
    INTEGER         ::     i,j

!!!$....................................................................
!!!$    
    !$OMP PARALLEL NUM_THREADS (ntd) DEFAULT(none) PRIVATE (cdum) &
    !$OMP SHARED (manz,nanz,sens,wmatd,wdfak,ap,bvec,lam,cgfac)
    !$OMP DO

    DO j=1,manz
       cdum = dcmplx(0d0)

       DO i=1,nanz
          cdum = cdum + dconjg(sens(i,j)) * &
               dcmplx(wmatd(i)*REAL(wdfak(i)))*ap(i)
       END DO

       bvec(j) = cdum + dcmplx(lam)*bvec(j)
       bvec(j) = bvec(j)*dcmplx(cgfac(j))
    END DO

    !$OMP END PARALLEL

  END SUBROUTINE bb


END MODULE cg_mod

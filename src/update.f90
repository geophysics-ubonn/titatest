SUBROUTINE update()
!!!$     
!!!$     Unterprogramm zum Bestimmen und Anbringen der Modellverbesserung
!!!$     mittels 'Smoothness Least Squares Method' und konjugierten
!!!$     Gradienten.
!!!$     Fuer beliebige Triangulierung und Stochastische Regularisierung
!!!$     
!!!$     Andreas Kemna                                            01-Mar-1996
!!!$     
!!!$     Letzte Aenderung                                         03-Aug-2009
!!!$     
!!!$.....................................................................

  USE alloci
  USE femmod
  USE datmod
  USE invmod
  USE cjgmod
  USE sigmamod
  USE modelmod
  USE elemmod
  USE cg_mod
  USE errmod
  USE konvmod

  IMPLICIT NONE


!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsvariablen
  COMPLEX (KIND(0D0)) ::   cdum
  REAL (KIND(0D0))    ::   dum,dum2

!!!$     Indexvariablen
  INTEGER (KIND=4)    ::  i,j,ij,in
!!!$.....................................................................

  IF (.NOT. llam) THEN

!!!$     Felder speichern
     dpar2 = dpar
     cgres2 = cgres

!!!$     Smoothnessvektor berechnen, bzw die RHS (Right Hand Side)
!!!$       b_q = A_q^h*C_d^-1*A_q*(d-f(m_q))-\lam C_m^-1 m_q
     IF (ltri==0) THEN
        DO i=1,manz
           cdum = dcmplx(0d0)
!!!$     diff+<
           IF (.NOT.lprior) THEN
!!!$ C_m^-1 * m_q
!!!$     diff+>
              IF (i.GT.1) cdum = dcmplx(smatm(i-1,2))*par(i-1)
              IF (i.LT.manz) cdum = cdum + dcmplx(smatm(i,2))*par(i+1)
              IF (i.GT.nx) cdum = cdum + dcmplx(smatm(i-nx,3))*par(i-nx)
              IF (i.LT.manz-nx+1) cdum = cdum + dcmplx(smatm(i,3))*par(i+nx)

              bvec(i) = cdum + dcmplx(smatm(i,1))*par(i)
!!!$     diff+<
           ELSE
!!!$ C_m^-1 * (m_q - m_0)
              IF (i.GT.1) cdum = dcmplx(smatm(i-1,2)) * (par(i-1)-m0(i-1))
              IF (i.LT.manz) cdum = cdum + dcmplx(smatm(i,2)) * &
                   (par(i+1)-m0(i+1))
              IF (i.GT.nx) cdum = cdum + dcmplx(smatm(i-nx,3)) * &
                   (par(i-nx)-m0(i-nx))
              IF (i.LT.manz-nx+1) cdum = cdum + dcmplx(smatm(i,3)) * &
                   (par(i+nx)-m0(i+nx))

              bvec(i)=cdum+dcmplx(smatm(i,1))*(par(i)-m0(i))

           END IF
!!!$     diff+>
        END DO
!!$c Damping--
     ELSE IF (ltri == 3.OR.ltri == 4) THEN

        WRITE(*,'(a)',ADVANCE='no')&
             'update:: damping has no part in the gradient'

!!!$     triang>
     ELSE IF (ltri == 1.OR.ltri == 2.OR. &
          (ltri > 4 .AND. ltri < 15)) THEN
!!!$ C_m^-1 * m_q
!!!$ C_m^-1 * (m_q - m_0)
        DO i=1,manz
           cdum = dcmplx(0d0)
           DO ij=1,smaxs
              in = nachbar(i,ij)
              IF (in /= 0) THEN
                 IF (.NOT. lprior) THEN
                    cdum = cdum + DCMPLX(smatm(i,ij)) * par(in)
                 ELSE
                    cdum = cdum + DCMPLX(smatm(i,ij)) * (par(in) - m0(in))
                 END IF
              END IF
           END DO
           IF (.NOT. lprior) THEN
              bvec(i) = cdum + DCMPLX(smatm(i,smaxs+1)) * par(i)
           ELSE
              bvec(i) = cdum + DCMPLX(smatm(i,smaxs+1)) * (par(i) - m0(i))
           END IF
        END DO

     ELSE IF (ltri == 15) THEN
        IF (.NOT. lprior) THEN

          !$OMP WORKSHARE
           bvec = MATMUL( DCMPLX(smatm),par)
           !$OMP END WORKSHARE

        ELSE

           !$OMP WORKSHARE
           bvec = MATMUL( DCMPLX(smatm),( par - m0 ) )
           !$OMP END WORKSHARE

        END IF
     END IF

!!!$ >> RM ref model regu
     IF (lw_ref) THEN
        DO i=1,manz
           IF (wref(i)) bvec(i) = bvec(i) + lam_ref * (par(i) - m_ref(i))
        END DO
     END IF
!!!$ << RM ref model regu

!!!$     triang<

!!!$  Skalierungsfaktoren bestimmen
!!!$ Preconditioning factors: 
!!!$ cgfac = diag(A^h_q * C_d^-1 * A + \lam C_m^-1)^-1
     DO j=1,manz
        dum = 0d0

        IF (ldc) THEN
           DO i=1,nanz
              dum = dum + sensdc(i,j) * sensdc(i,j) * wmatd(i) * &
                   DBLE(wdfak(i))
           END DO
        ELSE IF (lfpi) THEN
           DO i=1,nanz
              dum = dum + DBLE(sens(i,j)) * DBLE(sens(i,j)) * &
                   wmatd(i)*DBLE(wdfak(i))
           END DO
        ELSE
           DO i=1,nanz
              dum = dum + DBLE(dconjg(sens(i,j)) * sens(i,j)) * &
                   wmatd(i)*DBLE(wdfak(i))
           END DO
        END IF
        dum2 = dum

        IF (ltri==0) THEN

           dum    = dum + lam * smatm(j,1)

        ELSE IF (ltri == 1.OR.ltri == 2.OR. &
             (ltri > 4 .AND. ltri < 15)) THEN

           dum    = dum + lam * smatm(j,smaxs+1)

        ELSE IF (ltri == 3.OR.ltri == 4) THEN

           dum = dum + lam * smatm(j,1)

        ELSE IF (ltri == 15) THEN

           dum    = dum + lam * smatm(j,j)

        END IF

!!!$ >> RM ref model regu
        IF (wref(j)) dum = dum + lam * lam_ref
!!!!$<< RM

        cgfac(j) = 1d0/dsqrt(dum)
     END DO

!!!$     Konstantenvektor berechen und skalieren (RHS)
!!!$ the other part of the RHS system...
!!!$ A_q^h * C_d^-1 * (d-f(m_q)) + \lam C_m^-1 (m_q,(m_q-m_0))
    DO j=1,manz

        cdum = dcmplx(0d0)

!!!$     diff+<
        IF (.NOT.ldiff) THEN
!!!$     diff+>
           IF (ldc) THEN
              DO i=1,nanz
                 cdum = cdum + dcmplx(sensdc(i,j)*wmatd(i)* &
                      DBLE(wdfak(i)))*(dat(i)-sigmaa(i))
              END DO
           ELSE IF (lfpi) THEN
              DO i=1,nanz
                 cdum = cdum + dcmplx(DBLE(sens(i,j))*wmatd(i)* &
                      DBLE(wdfak(i)))*(dat(i)-sigmaa(i))
              END DO
           ELSE
              DO i=1,nanz
                 cdum = cdum + dconjg(sens(i,j))*dcmplx(wmatd(i)* &
                      DBLE(wdfak(i)))*(dat(i)-sigmaa(i))
              END DO
           END IF
!!!$     diff+<
        ELSE
           IF (ldc) THEN
              DO i=1,nanz
                 cdum = cdum + dcmplx(sensdc(i,j)*wmatd(i)* &
                      DBLE(wdfak(i)))*(dat(i)-sigmaa(i)-(d0(i)-fm0(i)))
              END DO
           ELSE IF (lfpi) THEN
              DO i=1,nanz
                 cdum = cdum + dcmplx(DBLE(sens(i,j))*wmatd(i)* &
                      DBLE(wdfak(i)))*(dat(i)-sigmaa(i)-(d0(i)-fm0(i)))
              END DO
           ELSE
              DO i=1,nanz
                 cdum = cdum + dconjg(sens(i,j))*dcmplx(wmatd(i)* &
                      DBLE(wdfak(i)))*(dat(i)-sigmaa(i)-(d0(i)-fm0(i)))
              END DO
           END IF
        END IF
!!!$     diff+>
        IF (ltri == 3.OR.ltri == 4) THEN

           bvec(j) = cdum

        ELSE
!!!$ RM ref model regu: reference model regu is already included in bvec (see above)..
           bvec(j) = cdum - dcmplx(lam)*bvec(j) 

        END IF

        bvec(j) = bvec(j)*dcmplx(cgfac(j))

!!$        print*,j,cdum,sens(1,j),sigma(j)

     END DO

!!$     DO i=1,nanz
!!$        print*,i,wmatd(i),dat(i),sigmaa(i),wdfak(i)
!!$     END Do
!!!$     Modellverbesserung mittels konjugierter Gradienten bestimmen
     CALL cjg

!!!$     Ggf. Verbesserung umspeichern
     IF (lfpi) THEN
        dpar = DCMPLX(0D0,DBLE(dpar))
     END IF

!!!$     Verbesserung skalieren
     dpar = dpar * DCMPLX(cgfac)
!!$     do j=1,manz
!!$        dpar(j) = dpar(j)*dcmplx(cgfac(j))
!!$     end do

  ELSE                      ! (llam==.TRUE.)


!!!$     Felder zuruecksetzen
     dpar = dpar2

     i = INT(cgres2(1))+1

     cgres(1:i) = cgres2(1:i)

  END IF

!!!$     Ggf. (Leitfaehigkeits-)Phasen < 0 mrad korrigieren
!!!$     if (lphi0.and.dimag(par(j)).lt.0d0)
!!!$     1        par(j) = dcmplx(dble(par(j)))
!!!$     akc Ggf. (Leitfaehigkeits-)Phasen < 1 mrad korrigieren
!!!$     ak            if (lphi0.and.dimag(par(j)).lt.1d-3)
!!!$     ak     1          par(j) = dcmplx(dble(par(j)),1d-3)

!!!$     i.e Stepsize = ||\delta m||
  bdpar = 0d0

  in = 0
  DO j=1,manz

     par(j) = par(j) + DCMPLX(step) * dpar(j) ! model update

!!$  ! eventually correct for phase < 0 mrad
     IF (lphi0 .AND. dimag(par(j)) < 0D0) THEN 
        par(j) = DCMPLX(DBLE(par(j))) 
        in = in + 1 
     END IF

     bdpar = bdpar + DBLE(dpar(j)*dconjg(dpar(j)))

  END DO

  IF (in > 0) WRITE (*,'(/a,I9,a/)',ADVANCE = 'no')' forcing zero ',in&
       ,' times'
  bdpar = bdpar * step

END SUBROUTINE update

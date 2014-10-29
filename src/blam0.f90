SUBROUTINE blam0()

!!!$     Unterprogramm zum Bestimmen des Start-Regularisierungsparameters.

!!!$     Andreas Kemna                                            20-Feb-1997
!!!$     Letzte Aenderung   07-Mar-2003

!!!$.....................................................................

  USE alloci
  USE femmod
  USE datmod
  USE invmod
  USE modelmod
  USE konvmod

  IMPLICIT NONE


!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsvariablen
  COMPLEX (KIND(0D0)) ::  cdum
  REAL (KIND(0D0))    ::  dum

  REAL (KIND(0D0)),ALLOCATABLE , DIMENSION(:) :: jtj

!!!$     Indexvariablen
  INTEGER (KIND = 4)  ::  i,j,k,ic

!!!$.....................................................................

!!!$     Start-Regularisierungsparameter bestimmen

!!!$ for fixed lambda set the values according to preset fixed lamfix
  IF (( BTEST(llamf,0) .OR. (lamnull_cri > EPSILON(lamnull_cri)) ) .AND..NOT. &
       lfpi ) THEN
     IF (nz==-1) THEN ! this is a special switch, but only taken for 
!!!!$ CRI/DC
        lammax = MAX(REAL(manz),REAL(nanz))
        WRITE (*,'(a,t5,a,G12.4)')ACHAR(13),'taking easy lam_0 ',lammax
     ELSE
        lammax = DBLE(lamnull_cri)
        WRITE (*,'(a,t5,a,G12.4)')ACHAR(13),'-> presetting lam0 CRI',lammax
     END IF
     RETURN
  ELSE IF ( BTEST(llamf,0) .OR. (lamnull_fpi > EPSILON(lamnull_fpi)) ) THEN
     lammax = DBLE(lamnull_fpi)
     WRITE (*,'(a,t5,a,G12.4)')ACHAR(13),'-> presetting lam0 FPI',lammax
     RETURN
  END IF


  ALLOCATE (jtj(manz))

  jtj = 0d0;ic = 0

  
  IF (ldc) THEN

     !$OMP PARALLEL DEFAULT(none) PRIVATE (dum) &
     !$OMP SHARED (manz,nanz,sensdc,wmatd,wdfak,jtj,lverb,ic)
     !$OMP DO SCHEDULE (GUIDED)

     DO j=1,manz
        IF (lverb) THEN
           !$OMP ATOMIC
           ic = ic + 1
           
           WRITE(*,'(a,t70,F6.2,A)',advance='no')ACHAR(13)//&
                'blam0/ ',REAL( ic * (100./manz)),'%'
        END IF
        dum = 0d0

        DO i=1,nanz
           DO k=1,manz
              dum = dum + sensdc(i,j) * sensdc(i,k) * &
                   wmatd(i)*DBLE(wdfak(i))
           END DO

        END DO

        jtj(j) = DABS(dum)
        
     END DO

     !$OMP END PARALLEL

  ELSE IF (lfpi) THEN

     !$OMP PARALLEL DEFAULT(none) PRIVATE (dum) &
     !$OMP SHARED (manz,nanz,sens,wmatd,wdfak,jtj,lverb,ic)
     !$OMP DO SCHEDULE (GUIDED)
     DO j=1,manz
        IF (lverb) THEN
           !$OMP ATOMIC
           ic = ic + 1
           
           WRITE(*,'(a,t70,F6.2,A)',advance='no')ACHAR(13)//&
                'blam0/ ',REAL( ic * (100./manz)),'%'
        END IF

        dum = 0d0

        DO i=1,nanz
           DO k=1,manz
              dum = dum + DBLE(sens(i,j)) * DBLE(sens(i,k)) * &
                   wmatd(i)*DBLE(wdfak(i))
           END DO
        END DO

        jtj(j) = DABS(dum)

     END DO
     !$OMP END PARALLEL

  ELSE

     !$OMP PARALLEL DEFAULT(none) PRIVATE (cdum) &
     !$OMP SHARED (manz,nanz,sens,wmatd,wdfak,jtj,lverb,ic)
     !$OMP DO SCHEDULE (GUIDED)
     DO j=1,manz

        
        IF (lverb) THEN
           !$OMP ATOMIC
           ic = ic + 1
           
           WRITE(*,'(a,t70,F6.2,A)',advance='no')ACHAR(13)//&
                'blam0/ ',REAL( ic * (100./manz)),'%'
        END IF

        cdum = dcmplx(0d0)

        DO i=1,nanz
           DO k=1,manz
              cdum = cdum + dconjg(sens(i,j)) * sens(i,k) * &
                   dcmplx(wmatd(i)*DBLE(wdfak(i)))
           END DO
        END DO
        
        jtj(j) = cdabs(cdum)

     END DO
     !$OMP END PARALLEL

  END IF

  lammax = SUM(jtj)/DBLE(manz)

  DEALLOCATE (jtj)

  lammax = lammax * 2d0/(alfx+alfz)
!!!$     ak Default
  lammax = lammax * 5d0
  WRITE (*,'(a,t5,a,G12.4,t60)')ACHAR(13),'found lam_0 ',lammax

!!!$     ak Synthetic Example (JoH)
!!!$     ak        lammax = lammax * 1d1

!!!$     ak MinFrac
!!!$     ak        lammax = lammax * 5d1

!!!$     ak Test
!!!$     ak        lammax = lammax * 1d1

!!!$     ak AAC
!!!$     ak        lammax = lammax * 5d0
  RETURN
END SUBROUTINE blam0

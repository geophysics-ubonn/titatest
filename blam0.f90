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

!!!$     Indexvariablen
  INTEGER (KIND = 4)  ::  i,j,k

!!!$.....................................................................

!!!$     Start-Regularisierungsparameter bestimmen

!!!$ for fixed lambda set the values according to preset fixed lamfix
  IF (( llamf .OR. (lamnull_cri > EPSILON(lamnull_cri)) ) .AND..NOT. &
       lip ) THEN
     IF (nz==-1) THEN ! this is a special switch, but only taken for 
!!!!$ CRI/DC
        lammax = MAX(REAL(manz),REAL(nanz))
        WRITE (*,'(t5,a,F12.1)')'taking easy lam_0 ',lammax
     ELSE
        lammax = DBLE(lamnull_cri)
        PRINT*,'-> presetting lam0 CRI',lammax
     END IF
     RETURN
  ELSE IF ( llamf .OR. (lamnull_fpi > EPSILON(lamnull_fpi)) ) THEN
     lammax = DBLE(lamnull_fpi)
     PRINT*,'-> presetting lam0 FPI',lammax
     RETURN
  END IF
  
  lammax = 0d0

  IF (ldc) THEN
     DO j=1,manz
        IF (lverb) WRITE(*,'(a,t70,F6.2,A)',advance='no')ACHAR(13)//&
             'blam0/ ',REAL( j * (100./manz)),'%'
        dum = 0d0

        DO i=1,nanz
           DO k=1,manz
              dum = dum + sensdc(i,j) * sensdc(i,k) * &
                   wmatd(i)*DBLE(wdfak(i))
           END DO
        END DO

        lammax = lammax + dabs(dum)
     END DO

  ELSE IF (lip) THEN

     DO j=1,manz
        IF (lverb) WRITE(*,'(a,t70,F6.2,A)',advance='no')ACHAR(13)//&
             'blam0/ ',REAL( j * (100./manz)),'%'
        dum = 0d0

        DO i=1,nanz
           DO k=1,manz
              dum = dum + DBLE(sens(i,j)) * DBLE(sens(i,k)) * &
                   wmatd(i)*DBLE(wdfak(i))
           END DO
        END DO

        lammax = lammax + dabs(dum)
     END DO

  ELSE

     DO j=1,manz
        IF (lverb) WRITE(*,'(a,t50,F6.2,A)',advance='no')ACHAR(13)//&
             'blam0/ ',REAL( j * (100./manz)),'%'
        cdum = dcmplx(0d0)

        DO i=1,nanz
           DO k=1,manz
              cdum = cdum + dconjg(sens(i,j)) * sens(i,k) * &
                   dcmplx(wmatd(i)*DBLE(wdfak(i)))
           END DO
        END DO

        lammax = lammax + cdabs(cdum)
     END DO

  END IF

  lammax = lammax/DBLE(manz)

  lammax = lammax * 2d0/(alfx+alfz)
!!!$     ak Default
  lammax = lammax * 5d0
  IF (lammax > 1e10) THEN
     WRITE (*,'(t5,a,G12.1)')'found lam_0 ',lammax
  ELSE
     WRITE (*,'(t5,a,F12.1)')'found lam_0 ',lammax
  END IF
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

subroutine blam0()

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

  IMPLICIT none


!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsvariablen
  COMPLEX (KIND(0D0)) ::  cdum
  REAL (KIND(0D0))    ::   dum,lam_tmp

!!!$     Indexvariablen
  INTEGER (KIND = 4)  ::  i,j,k,count

!!!$.....................................................................

!!!$     Start-Regularisierungsparameter bestimmen
  IF (nz < 0.OR.lamfix /= 0d0) THEN
     lammax = ABS(lamfix)
     IF (nz==-1) lammax = MAX(REAL(manz),REAL(nanz))
     WRITE (*,'(t5,a,F12.1)')'taking easy lam_0 ',lammax
     RETURN
  END IF

  lammax = 0d0
  count = 0
  !$OMP PARALLEL DEFAULT (none) &
  !$OMP SHARED (manz,count,lverb,ldc,nanz,sensdc,wmatd,wdfak,lip,sens,lammax) &
  !$OMP PRIVATE(i,k,dum,j,cdum)
  !$OMP DO
  do j=1,manz

     !$OMP ATOMIC
     count = count + 1

     IF (lverb) write(*,'(a,t70,F6.2,A)',advance='no')ACHAR(13)//&
          'blam0/ ',REAL( count * (100./manz)),'%'
     dum = 0d0;cdum = DCMPLX(0D0)
     if (ldc) then
        do i=1,nanz
           do k=1,manz
              dum = dum + sensdc(i,j) * sensdc(i,k) * &
                   wmatd(i)*dble(wdfak(i))
           end do
        end do
     else if (lip) then
        do i=1,nanz
           do k=1,manz
              dum = dum + dble(sens(i,j)) * dble(sens(i,k)) * &
                   wmatd(i)*dble(wdfak(i))
           end do
        end do
     else
        do i=1,nanz
           do k=1,manz
              cdum = cdum + dconjg(sens(i,j)) * sens(i,k) * &
                   dcmplx(wmatd(i)*dble(wdfak(i)))
           end do
        end do
        dum = CDABS(cdum)
     END if
     
     !$OMP ATOMIC
     lammax = lammax + dabs(dum)

  end do
  !$OMP END PARALLEL

  lammax = lammax/dble(manz)

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
end subroutine blam0

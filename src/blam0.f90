subroutine blam0()

!!!$     Unterprogramm zum Bestimmen des Start-Regularisierungsparameters.

!!!$     Andreas Kemna                                            20-Feb-1997
!!!$     Letzte Aenderung   07-Mar-2003

!!!$.....................................................................

  use alloci
  use femmod
  use datmod
  use invmod
  use modelmod
  use konvmod

  implicit none


!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsvariablen
  complex (prec) ::  cdum,zlange
  real (prec)    ::  dum
  complex (prec),dimension(:,:),allocatable :: WDA,awdwda,wdah
  complex(prec),dimension(:,:),allocatable:: tmp
  real (prec),allocatable , dimension(:) :: jtj

!!!$     Indexvariablen
  integer (KIND = 4)  ::  i,j,k,ic

!!!$.....................................................................

!!!$     Start-Regularisierungsparameter bestimmen

!!!$ for fixed lambda set the values according to preset fixed lamfix
  if (( btest(llamf,0) .or. (lamnull_cri > epsilon(lamnull_cri)) ) .and..not. &
       lfpi ) then
     if (nz==-1) then ! this is a special switch, but only taken for 
!!!!$ CRI/DC
        lammax = max(real(manz),real(nanz))
        write (*,'(a,t5,a,G12.4)')achar(13),'taking easy lam_0 ',lammax
     else
        lammax = real(lamnull_cri)
        write (*,'(a,t5,a,G12.4)')achar(13),'-> presetting lam0 CRI',lammax
     end if
     return
  else if ( btest(llamf,0) .or. (lamnull_fpi > epsilon(lamnull_fpi)) ) then
     lammax = real(lamnull_fpi)
     write (*,'(a,t5,a,G12.4)')achar(13),'-> presetting lam0 FPI',lammax
     return
  end if


  allocate (jtj(manz))

  jtj = 0d0;ic = 0


  if (ldc) then

     !$OMP PARALLEL DEFAULT(none) PRIVATE (dum) &
     !$OMP SHARED (manz,nanz,sensdc,wmatd,wdfak,jtj,lverb,ic)
     !$OMP DO SCHEDULE (GUIDED)

     do j=1,manz
        if (lverb) then
           !$OMP ATOMIC
           ic = ic + 1

           write(*,'(a,t70,F6.2,A)',advance='no')achar(13)//&
                'blam0/ ',real( ic * (100./manz)),'%'
        end if
        dum = 0d0

        do i=1,nanz
           do k=1,manz
              dum = dum + sensdc(i,j) * sensdc(i,k) * &
                   wmatd(i)*real(wdfak(i))
           end do

        end do

        jtj(j) = abs(dum)

     end do

     !$OMP END PARALLEL

  else if (lfpi) then

     !$OMP PARALLEL DEFAULT(none) PRIVATE (dum) &
     !$OMP SHARED (manz,nanz,sens,wmatd,wdfak,jtj,lverb,ic)
     !$OMP DO SCHEDULE (GUIDED)
     do j=1,manz
        if (lverb) then
           !$OMP ATOMIC
           ic = ic + 1

           write(*,'(a,t70,F6.2,A)',advance='no')achar(13)//&
                'blam0/ ',real( ic * (100./manz)),'%'
        end if

        dum = 0d0

        do i=1,nanz
           do k=1,manz
              dum = dum + real(sens(i,j)) * real(sens(i,k)) * &
                   wmatd(i)*real(wdfak(i))
           end do
        end do

        jtj(j) = abs(dum)

     end do
     !$OMP END PARALLEL

  else
! Perform lambda_0 search. See Kemna (2000), page 67.
! complex case  
  
     allocate(wda(nanz,manz), wdah(manz,nanz))
     allocate(awdwda(manz,manz))
     lammax = 0.
     wdah = conjg( transpose( sens ) )
     do i = 1,nanz
        wda(i,:) = cmplx(wmatd(i)*real(wdfak(i)))*sens(i,:)
     end do
     call zgemm('n', 'n', manz, manz, nanz, 1D0, wdah, manz, wda, nanz, 0D0,&
        awdwda,manz)
     jtj  = abs(sum((awdwda),1))
     lammax = sum((jtj),1)
     
     deallocate(wda,wdah,awdwda)
  end if
  
  lammax = lammax/real(manz)

  deallocate (jtj)

  lammax = lammax * 2d0/(alfx+alfz)
!!!$     ak Default
  lammax = lammax * 5d0
  write (*,'(F13.2)') lammax

!!!$     ak Synthetic Example (JoH)
!!!$     ak        lammax = lammax * 1d1

!!!$     ak MinFrac
!!!$     ak        lammax = lammax * 5d1

!!!$     ak Test
!!!$     ak        lammax = lammax * 1d1

!!!$     ak AAC
!!!$     ak        lammax = lammax * 5d0
  return
end subroutine blam0

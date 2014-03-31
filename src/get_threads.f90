SUBROUTINE GET_THREADS(nthreads,kwnanz)

!Determine the optimum number of threads for CRTomo. 
!It's NOT always the CPU count :)

USE omp_lib
IMPLICIT NONE
! Output
integer :: nthreads

! Input
integer :: kwnanz

! Internal
integer :: maxthreads, mythreads, i
############################################
maxthreads = 1
!maxthreads = OMP_GET_MAX_THREADS()
NTHREADS = maxthreads
! now that we know nf and kwnanz, we can adjust the OMP environment..
IF (maxthreads > 2) THEN 
  ! single or double processor machines don't need scheduling
  mythreads = MAX(kwnanz,2)
  IF ( mythreads <= maxthreads ) THEN ! best case,
    ! the number of processors is greater or equal the assumed workload
  ELSE 
    ! is smaller than the minimum workload. now we have to devide a bit
    DO i = 1, INT(kwnanz/2)
      mythreads = INT(kwnanz / i) + 1
      IF (mythreads < maxthreads) EXIT
    END DO
  END IF
  NTHREADS = mythreads
END IF
!CALL OMP_SET_NUM_THREADS ( NTHREADS )
! recheck ..
!i = OMP_GET_MAX_THREADS()
WRITE(6,'(2(a, i3),a)') " openMP threads: ",i,'(',maxthreads,' CPUs)'
END SUBROUTINE GET_THREADS

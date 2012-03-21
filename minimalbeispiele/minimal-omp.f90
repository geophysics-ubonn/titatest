MODULE my_vars


  REAL,PUBLIC,DIMENSION (:),ALLOCATABLE :: z
  REAL,PUBLIC :: r
  REAL,PUBLIC :: pi,norm_r,norm_z
  REAL,PUBLIC :: fill

  INTEGER,PUBLIC :: i
  INTEGER,PUBLIC  :: TID ! thread ID
  INTEGER,PUBLIC  :: NTHREADS
  INTEGER,PUBLIC  :: n

END MODULE my_vars

PROGRAM minimal_parallel
  USE omp_lib
  USE my_vars

  IMPLICIT none


  pi = ACOS(-1.)

  n=AINT(pi)
  n=1
  print *,n,pi
  write(6,"(a, i3)") " OpenMP max threads: ", OMP_GET_MAX_THREADS()

  ALLOCATE (z(OMP_GET_MAX_THREADS()))
!  ALLOCATE (r(n))

  !$OMP PARALLEL PRIVATE(TID,r) SHARED(z)
  r = 1.0
  TID = OMP_GET_THREAD_NUM()
  PRINT*,'Hello from thread =', TID
!!$ this is for the master thread
  IF (TID == 0) THEN
     NTHREADS = OMP_GET_NUM_THREADS()
     write(6,'(a,i3)') " OpenMP master: N_threads = ",NTHREADS
     r = 1.0
  ELSE
     r = pi
  END IF
  print*,'main call :: my_sub::',omp_get_thread_num(),r,z

  CALL my_sub(r)!,z,OMP_GET_MAX_THREADS())

!!!$ all threads rejoin master thread and disband
  !$OMP END PARALLEL
  
  print*,'main program z',z
  
END PROGRAM minimal_parallel

SUBROUTINE my_sub(my_r)!,z,n)

  USE omp_lib
  USE my_vars
  IMPLICIT none

  REAL,INTENT(INOUT) :: my_r

!!$  INTEGER :: n
!!$  REAL,INTENT(INOUT),DIMENSION(n) :: z

  PRINT*,'on_entry::',omp_get_thread_num(),my_r

  my_r = my_r + 1
!  r = r + 1
  z(omp_get_thread_num()+1) = z(omp_get_thread_num()+1) + my_r

  PRINT*,'on_return::',omp_get_thread_num(),my_r

END SUBROUTINE my_sub

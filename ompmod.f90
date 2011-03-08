MODULE ompmod

  IMPLICIT none

  INTEGER,PUBLIC :: TID ! Thread ID
  INTEGER,PUBLIC :: NTHREADS ! Total number of threads
  INTEGER,PUBLIC :: CHUNK ! chunk pieces for static scheduling of do loops, a.k.
  INTEGER,EXTERNAL :: OMP_GET_MAX_THREADS !OMP function
  INTEGER,EXTERNAL :: OMP_GET_NUM_THREADS ! OMP function
  INTEGER,EXTERNAL :: OMP_GET_THREAD_NUM ! OMP function

  
!!$  !$OMP PARALLEL SHARED(a,b,CHUNK) PRIVATE (i)
!!$  a = 1.
!!$  !$OMP DO SCHEDULE(STATIC, chunk)

!!$  DO i = 1, n
!!$    a(i) = a(i) + b(i)
!!$  END DO

!!$  !$OMP END PARALLEL

END MODULE ompmod

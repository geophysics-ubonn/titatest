SUBROUTINE chold(a,p,n,ierr,lverb)
!!$c-----------------------------------------------------------------
!!$c
!!$c                      Cholesky Decomposition
!!$c                      **********************
!!$c
!!$c This subroutine calculates the lower triangular matrix L which,
!!$c when multiplied by its own transpose, gives the symmetric 
!!$c matrix A. From Numerical recipies (Press et al 2003)
!!$c 
!!$c Changed and put into this format by R. Martin 2010
!!$c
!!$c**********          NOTE           ***************************
!!$c chold In this form needs only the upper triangular part of a
!!$c and stores L in the lower part.
!!$c P contains the diagonal entries (EVs)
!!$c**************************************************************
!!$c INPUT VARIABLES:
!!$c
!!$c   a(n,n) Symmetric positive definite nxn matrix
!!$c          to be decomposed
!!$c          Upper part still contains A and lower part is filled with L
!!$c   p(n)   Eigenvalues of a 
!!$c   ierr   Error code:  ierr=0 - no errors; ierr=1 - matrix a
!!$c          is not positive definite
!!$c NO EXTERNAL REFERENCES:
!!$c------------------------------------------------------------
  USE ompmod

  IMPLICIT none

  INTEGER,INTENT (IN)                             :: n
  REAL (KIND(0D0)), DIMENSION (n,n),INTENT(INOUT) :: a
  REAL (KIND(0D0)), DIMENSION (n),INTENT(OUT)     :: p
  LOGICAL,INTENT(IN)                              :: lverb
  INTEGER, INTENT (OUT)                           :: ierr
  REAL (KIND(0D0))                                :: s
  INTEGER                                         :: i,k,j,count

  ierr = 0
  count = 0
!!$  !$OMP PARALLEL DEFAULT (none) &
!!$  !$OMP PRIVATE (j,s,k) &
!!$  !$OMP SHARED (a,p,n,ierr,count)
!!$  !$OMP DO ORDERED
  DO i = 1 , n
     
     !$OMP ATOMIC
     count = count + 1
     
     IF (lverb) WRITE (*,'(A,t25,F6.2,A)',ADVANCE='no')&
          ACHAR(13)//'Factorization',REAL( count * (100./n)),'%'
     
     DO j = i , n
        
        s = a(i,j)
        
        DO k = i-1 , 1 ,-1
           
           s = s - a(i,k) * a(j,k) ! line sum
           
        END DO
        
        IF (i == j) THEN
           
           IF (s <= 0) THEN
              PRINT*,'CHOLD:: - not positive definite', s
              ierr = -i
              STOP
           END IF

           p(i) = DSQRT(s) ! main diagonal

        ELSE

           a(j,i) = s / p(i) ! scale value

        END IF
        
     END DO
  END DO

!!$ !$OMP END PARALLEL

END SUBROUTINE chold

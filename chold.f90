SUBROUTINE chold(a,n,u,ierr)
!!$c-----------------------------------------------------------------
!!$c
!!$c                      Cholesky Decomposition
!!$c                      **********************
!!$c
!!$c This subroutine calculates the lower triangular matrix L which,
!!$c when multiplied by its own transpose, gives the symmetric 
!!$c matrix A. From Numerical recipies
!!$c 
!!$c**********          NOTE           ***************************
!!$c chold In this form needs only the upper triangular part of a
!!$c and stores L in the lower part.
!!$c P contains the diagonal entries (EVs)
!!$c**************************************************************
!!$c INPUT VARIABLES:
!!$c
!!$c   a(n,n) Symmetric positive definite nxn matrix
!!$c    to be decomposed 
!!$c   u(n,n) upper triangle with U^TU = A
!!$c   ierr   Error code:  ierr=0 - no errors; ierr=1 - matrix a
!!$c          is not positive definite
!!$c NO EXTERNAL REFERENCES:
!!$c------------------------------------------------------------
  IMPLICIT none

  INTEGER,INTENT (IN)               :: n
  REAL (KIND(0D0)), DIMENSION (n,n) :: a,u
  REAL (KIND(0D0))                  :: s,d
  REAL (KIND(0D0)),PARAMETER        :: tol=1.e-10
  INTEGER, INTENT (OUT)             :: ierr
  INTEGER                           :: i,k,j

  ierr = 0
  u = 0.

  DO i = 1,n
     WRITE (*,'(A,1X,F6.2,A)',ADVANCE='no')&
          ACHAR(13)//ACHAR(9)//ACHAR(9)//&
          ACHAR(9)//'/ ',REAL( i * (100./n)),'%'
     s = DOT_PRODUCT(u(1:i,i),u(1:i,i))
     d = a(i,i) - s
     IF (ABS(d)< tol) THEN
        u(i,i) = 0.
     ELSE
        IF (d < 0.) THEN
           PRINT*,'WARNING: chold - not positive definite'
           ierr = -i
           RETURN
        END IF
        u(i,i) = SQRT(d)
     END IF

     DO j = i+1,n
        s = DOT_PRODUCT(u(1:i,i),u(1:i,j))
        IF (ABS(s) < tol) s = 0.
        u(i,j) = (a(i,j) - s) / u(i,i)
     END DO

  END DO

  ierr = 0

END SUBROUTINE chold

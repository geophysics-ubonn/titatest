SUBROUTINE chold(a,p,n,ierr)
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
  IMPLICIT none

  INTEGER,INTENT (IN)               :: n
  REAL (KIND(0D0)), DIMENSION (n,n) :: a
  REAL (KIND(0D0)), DIMENSION (n)   :: p
  REAL (KIND(0D0))                  :: s
  INTEGER, INTENT (OUT)             :: ierr
  INTEGER                           :: i,k,j

  ierr = 0

  DO i = 1 , n

     WRITE (*,'(A,1X,F6.2,A)',ADVANCE='no')&
          ACHAR(13)//ACHAR(9)//ACHAR(9)//&
          ACHAR(9)//'/ ',REAL( i * (100./n)),'%'

     DO j = i , n

        s = a(i,j)

        DO k = i-1 , 1 ,-1

           s = s - a(i,k) * a(j,k) ! line sum

        END DO

        IF (i == j) THEN

           IF (s <= 0) THEN
              PRINT*,'CHOLD:: - not positive definite', s
              ierr = -i
              RETURN
           END IF

           p(i) = DSQRT(s) ! main diagonal

        ELSE

           a(j,i) = s / p(i) ! scale value

        END IF

     END DO
  END DO

  ierr = 0

END SUBROUTINE chold

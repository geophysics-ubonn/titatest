SUBROUTINE linv(a,b,n)
!!$c-----------------------------------------------------------------------
!!$c
!!$c                Inverse of a Lower Triangular Matrix
!!$c                ************************************
!!$c
!!$c This subroutine finds the inverse of a lower triangular matrix A and
!!$c stores the answer in B. (from "Numerical Analysis of Symmetric
!!$c Matrices,"  H.R. Schwarz et al.,)
!!$c
!!$c
!!$c
!!$c INPUT VARIABLES:
!!$c
!!$c   a(n,n)           Lower triangular matrix to be inverted
!!$c   b(n,n)           the inverse
!!$c   n                Dimension of the matrix you're inverting
!!$c-----------------------------------------------------------------------
  IMPLICIT none

  INTEGER,INTENT (IN)               :: n
  REAL (KIND(0D0)), DIMENSION (n,n) :: a,b
  REAL (KIND(0D0))                  :: sum
  INTEGER                           :: k,i,j
  
  b = 0.
  DO i = 1,n
     IF(i>1) THEN
        DO k = 1,i-1
           sum=0.
           DO j = k,i-1
              sum = sum + a(i,j)*b(j,k)
           END DO
           b(i,k) = -sum/a(i,i)
        END DO
     END IF
     b(i,i) = 1./a(i,i)
  END DO
!!$c
!!$c Finished:
!!$c

END SUBROUTINE linv
    

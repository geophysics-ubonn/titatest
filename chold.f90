SUBROUTINE chold(a,t,n,ierr)

!!$c-----------------------------------------------------------------------
!!$c
!!$c                      Cholesky Decomposition
!!$c                      **********************
!!$c
!!$c This subroutine calculates the lower triangular matrix T which, when
!!$c multiplied by its own transpose, gives the symmetric matrix A. (from
!!$c "Numerical Analysis of Symmetric Matrices,"  H.R. Schwarz et al.,
!!$c p. 254)
!!$c
!!$c
!!$c
!!$c INPUT VARIABLES:
!!$c
!!$c   a(n,n)           Symmetric positive definite matrix to be
!!$c                      decomposed (destroyed in the calculation of t)
!!$c   t(n,n)           Lower triangular matrix solution
!!$c   n                Dimension of the system you're decomposing
!!$c   ierr             Error code:  ierr=0 - no errors; ierr=1 - matrix a
!!$c                      is not positive definite
!!$c
!!$c
!!$c
!!$c NO EXTERNAL REFERENCES:
!!$c-----------------------------------------------------------------------
  IMPLICIT none

  INTEGER,INTENT (IN)               :: n
  REAL (KIND(0D0)), DIMENSION (n,n) :: a,t
  INTEGER, INTENT (OUT)             :: ierr
  INTEGER                           :: ip,k,i

  ierr = 0
!!$c
!!$c Check for positive definiteness:
!!$c
  DO ip=1,n
     IF(a(ip,ip)<=0.0) THEN
        PRINT*,'WARNING: chol - not positive definite'
        ierr = 1
        RETURN
     END IF
     t(ip,ip) = SQRT (a(ip,ip))
     IF(ip>=n) RETURN
     DO k = ip+1,n
        t(k,ip) = a(ip,k)/t(ip,ip)
     END DO
     DO i = ip+1,n
        DO k = i,n
           a(i,k) = a(i,k) - t(i,ip) * t(k,ip)
        END DO
     END DO
  END DO
!!$c
!!$c Finished
!!$c
END SUBROUTINE chold

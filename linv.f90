SUBROUTINE linv(u,b,n)
!!$c-----------------------------------------------------------------
!!$c
!!$c                Inverse of a Lower Triangular Matrix
!!$c                ************************************
!!$c
!!$c This subroutine finds the inverse of factorized Lower triangular
!!$c 
!!$c from Numerical recipies
!!$c
!!$c
!!!$ in order to find the inverse of A we need to compute 
!!!$ (U^-1)^T U^-1 , since
!!!$ (U^-1)^T U^-1 U U^T = I
!!!$
!!$c INPUT VARIABLES:
!!$c
!!$c   u(n,n)  Upper triangular matrix from chold
!!$c   b(n,n)  Inverse of former A
!!$c   n       leading dimension 
!!$c-----------------------------------------------------------------
  IMPLICIT none
  
  INTEGER,INTENT (IN)               :: n
  REAL (KIND(0D0)), DIMENSION (n,n) :: u,b
  REAL (KIND(0D0))                  :: s,ujj
  INTEGER                           :: k,i,j
  
  b = 0.

  DO j = n,1,-1
     WRITE (*,'(A,1X,F6.2,A)',ADVANCE='no')&
          ACHAR(13)//ACHAR(9)//ACHAR(9)//&
          ACHAR(9)//'/ ',REAL( (n-j) * (100./n)),'%'
     ujj = u(j,j)
     s = DOT_PRODUCT(u(j,j+1:n),b(j,j+1:n))
     b(j,j) = 1. / ujj**2. - s / ujj
     DO i = j,1,-1
        s = DOT_PRODUCT(u(i,i+1:n),b(i+1:n,j))
        b(i,j) = - s / u(i,i)
        b(j,i) = b(i,j)
     END DO
  END DO

END SUBROUTINE linv
    

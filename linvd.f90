SUBROUTINE linvd(a,p,n)
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
!!!$ in order to find the inverse of A we still need to compute 
!!!$ (L^-1)^T L^-1 , since
!!!$ (L^-1)^T L^-1 L L^T = I
!!!$
!!$c INPUT VARIABLES:
!!$c
!!$c   a(n,n)  Lower triangular matrix from chold 
!!$c           which stores its inverse of a (!) in the upper part..
!!$c   p(n)    Eigenvalues of former A
!!$c   n       leading dimension 
!!$c-----------------------------------------------------------------
  IMPLICIT none
  
  INTEGER,INTENT (IN)               :: n
  REAL (KIND(0D0)), DIMENSION (n,n) :: a
  REAL (KIND(0D0)), DIMENSION (n)   :: p
  REAL (KIND(0D0))                  :: s
  INTEGER                           :: k,i,j
  
  DO i= 1 , n
     WRITE (*,'(A,1X,F6.2,A)',ADVANCE='no')&
          ACHAR(13)//ACHAR(9)//ACHAR(9)//&
          ACHAR(9)//'/ ',REAL( (n-i) * (100./n)),'%'

     a(i,i) = 1. / p(i)

     DO j = i+1 , n

        s = 0.

        DO k = i , j-1

           s = s - a(j,k) * a(k,i)

        END DO

        a(j,i) = s / p(j)
        
        a(i,j) = 0. ! upper triangle eq zero

     END DO

  END DO

  DO i = 1, n


     a(i,i) = a(i,i) * a(i,i) 

     DO k = i+1 , n

        a(i,i) = a(i,i) + a(k,i) * a(k,i)

     END DO
     
     DO j = i + 1, n
        
        DO k = j, n

           a(i,j) = a(i,j) + a(k,i) * a(k,j)

        END DO
     END DO

  END DO

END SUBROUTINE linvd
    

SUBROUTINE Gauss (a,n,e_flag)    ! Invert matrix by Gauss method
  IMPLICIT NONE
  INTEGER,INTENT (IN)            :: n
  REAL(KIND(0D0)),DIMENSION(n,n) :: a,b
  REAL(KIND(0D0)),DIMENSION(n)   :: temp
  REAL(KIND(0D0))                :: c,d 
  INTEGER                        :: j,k,m,imax(1),e_flag
  INTEGER,DIMENSION (n)          :: ipvt
  
  b = a
  
  e_flag=-1
  
  DO j=1,n
     ipvt(j) = j
  END DO
  
  DO k = 1,n
     imax = MAXLOC(ABS(b(k:n,k)))
     m = k-1+imax(1)
     IF (m /= k) THEN
        PRINT*,'Pivoting:: ',ipvt( [m,k] )
        ipvt( [m,k] ) = ipvt( [k,m] )
        b( [m,k],:) = b( [k,m],:)
     END IF
     d = 1/b(k,k)
     temp = b(:,k)
     DO j = 1, n
        c = b(k,j)*d
        b(:,j) = b(:,j)-temp*c
        b(k,j) = c
     END DO
     b(:,k) = temp*(-d)
     b(k,k) = d
  END DO
  a(:,ipvt) = b
  
  e_flag=0
  
END SUBROUTINE Gauss



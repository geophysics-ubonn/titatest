!     --------------------------------------------------------------------
      SUBROUTINE Gauss (a,n,e_flag)    ! Invert matrix by Gauss method
!     --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: n
      REAL(8) :: a(n,n), b(n,n), c, d, temp(n)
      INTEGER :: j, k, m, imax(1), ipvt(n),e_flag

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


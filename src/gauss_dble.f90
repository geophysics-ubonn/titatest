SUBROUTINE Gauss_REAL (a,n,e_flag)    ! Invert matrix by Gauss method
use alloci, only: prec
  IMPLICIT NONE
  INTEGER(KIND(4))                             :: n
  REAL(prec),DIMENSION(n,n),INTENT(INOUT) :: a
  REAL(prec),DIMENSION(:), ALLOCATABLE    :: temp
  REAL(prec)                              :: c,d 
  INTEGER(KIND(4))                             :: i,j,e_flag
  

  ALLOCATE (temp(n),STAT=e_flag)
  IF (e_flag/=0) THEN
     print*,'error alllocating temp(',n,')=',n*8/(1024**3),' GB'
     RETURN
  END IF

  DO i = 1,n
     WRITE (*,'(A,1X,F6.2,A)',ADVANCE='no')ACHAR(13)//ACHAR(9)//ACHAR(9)//&
     ACHAR(9)//'/ ',REAL( i * (100./n)),'%'

     e_flag = -i

     IF (ABS(a(i,i)) < EPSILON(d)) RETURN

     d = 1D0/a(i,i)
     temp = a(:,i)
     DO j = 1, n
        c = a(i,j)*d
        a(:,j) = a(:,j)-temp*c
        a(i,j) = c
     END DO
     a(:,i) = temp*(-d)
     a(i,i) = d
  END DO

  e_flag=0
  
END SUBROUTINE Gauss_REAL

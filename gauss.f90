SUBROUTINE Gauss (n,e_flag)    ! Invert matrix by Gauss method

  USE alloci,only:smatm

  IMPLICIT NONE
  INTEGER(KIND(4))                             :: n
  REAL(KIND(0D0)),DIMENSION(:,:), ALLOCATABLE  :: dump
  REAL(KIND(0D0)),DIMENSION(:), ALLOCATABLE    :: temp
  REAL(KIND(0D0))                              :: c,d 
  INTEGER                                      :: j,k,m,imax(1),e_flag
  INTEGER,DIMENSION (:), ALLOCATABLE           :: ipvt
  

  ALLOCATE (dump(n,n),STAT=e_flag)
  
  IF (e_flag/=0) THEN
     print*,'error alllocating dump(',n,',',n,')=',n*n*16/(1024**3),' GB'
     STOP
  END IF
  
  dump = smatm

  ALLOCATE (temp(n),ipvt(n),STAT=e_flag)

  e_flag=-1
  
  DO j=1,n
     ipvt(j) = j
  END DO
  
  DO k = 1,n
     WRITE (*,'(A,1X,F5.2,A)',ADVANCE='no')ACHAR(13)//'gauss/ ',&
          REAL( k * (100/n)),'%'
     imax = MAXLOC(ABS(dump(k:n,k)))
     m = k-1+imax(1)

     IF (m /= k) THEN
        PRINT*,'Pivoting:: ',ipvt( [m,k] )
        ipvt( [m,k] ) = ipvt( [k,m] )
        dump( [m,k],:) = dump( [k,m],:)
     END IF

     d = 1/dump(k,k)
     temp = dump(:,k)
     DO j = 1, n
        c = dump(k,j)*d
        dump(:,j) = dump(:,j)-temp*c
        dump(k,j) = c
     END DO
     dump(:,k) = temp*(-d)
     dump(k,k) = d

  END DO

  smatm = dump

  DEALLOCATE (dump,temp,ipvt)

  e_flag=0
  
END SUBROUTINE Gauss



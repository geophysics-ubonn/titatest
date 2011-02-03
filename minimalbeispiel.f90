PROGRAM minimal
  IMPLICIT none
  COMPLEX(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: dz_vct
  COMPLEX(KIND(0D0))                            :: dz
  REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: r_vct
  REAL(KIND(0D0))                            :: r,pi
  INTEGER                              :: i,n


  n = 30

  pi = ACOS(-1d0)
  print*,'pi = ACOS(-1d0):',pi
  pi = ACOS(-1.)
  print*,'pi = ACOS(-1.):',pi
  pi = DACOS(-1d0)
  print*,'pi = DACOS(-1d0):',pi

  ALLOCATE (dz_vct(n))


  
END PROGRAM minimal

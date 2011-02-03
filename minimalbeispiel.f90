PROGRAM main
  IMPLICIT none
  COMPLEX(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: dz_vct
  COMPLEX(KIND(0D0))                            :: dz,dz2
  REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: r_vct
  REAL(KIND(0D0))                            :: r,pi
  INTEGER                              :: i,n


  n = 100
  print*,'KIND(pi):',KIND(pi)
  print*
  pi = ACOS(-1d0)
  print*,'pi = ACOS(-1d0):',pi
  pi = ACOS(-1.)
  print*,'pi = ACOS(-1.):',pi
  pi = DACOS(-1d0)
  print*,'pi = DACOS(-1d0):',pi
  r = SQRT(pi)
  print*,'r = SQRT(pi):',r
  r = DSQRT(pi)
  print*,'r = DSQRT(pi):',r
  print*
  print*,'KIND(dz):',KIND(dz)
  print*
  dz = CMPLX(1.0,pi)
  print*,'dz = CMPLX(1.0,pi):',dz
  dz = CMPLX(pi,pi)
  print*,'dz = CMPLX(pi,pi):',dz
  dz = DCMPLX(pi,pi)
  print*,'dz = DCMPLX(pi,pi):',dz
  dz = DCMPLX(1.0,pi)
  print*,'dz = DCMPLX(1.0,pi):',dz
  print*
  dz = DCMPLX(pi,pi)
  print*,'dz = DCMPLX(pi,pi):',dz
  dz2 = SQRT(dz)
  print*,'dz2 = SQRT(dz):',dz2
  dz2 = CDSQRT(dz)
  print*,'dz2 = CDSQRT(dz):',dz2
  dz2 = EXP(dz)
  print*,'dz2 = EXP(dz):',dz2
  dz2 = CDEXP(dz)
  print*,'dz2 = CDEXP(dz):',dz2
  dz2 = EXP(DCMPLX(pi,pi))
  print*,'dz = EXP(DCMPLX(pi,pi)):',dz2
  dz2 = CDEXP(DCMPLX(pi,pi))
  print*,'dz = CDEXP(DCMPLX(pi,pi)):',dz2
  print*

  ALLOCATE (dz_vct(n),r_vct(n))

  r_vct = pi
  print*,'KIND(r_vct):',KIND(r_vct)
  print*,'SHAPE(r_vct):',SHAPE(r_vct)
  DO i=1,n/3,2
     print*,i,r_vct(i)
  END DO

  r=0d0
  DO i=1,n
     r = r + r_vct(i)*r_vct(i)
  END DO
  print*,'Inner product of r_vct via do loop',r
  print*,'and via DOT_PRODUCT(r_vct,r_vct)',DOT_PRODUCT(r_vct,r_vct)
  print*

  print*,'KIND(dz_vct):',KIND(dz_vct)
  print*,'SHAPE(dz_vct):',SHAPE(dz_vct)

  dz_vct = DCMPLX(pi,pi)

  DO i=1,n/3,2
     print*,i,dz_vct(i)
  END DO
  r=0d0
  DO i=1,n
     r = r + (CONJG(dz_vct(i))*dz_vct(i))
  END DO
  print*,'Inner product of dz_vct via do loop',r
  print*,'and via DOT_PRODUCT(CONJG(dz_vct),dz_vct)',&
       DOT_PRODUCT(CONJG(dz_vct),dz_vct)
  print*,'AIMAG(DOT_PRODUCT(CONJG(dz_vct),dz_vct))',&
       AIMAG(DOT_PRODUCT(CONJG(dz_vct),dz_vct))
  print*
  
END PROGRAM main

PROGRAM main
  IMPLICIT none
  COMPLEX(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: dz_vct
  COMPLEX(KIND(0D0))                            :: dz,dz2
  COMPLEX(KIND(0E0)),DIMENSION(:),ALLOCATABLE   :: z_vct
  COMPLEX(KIND(0E0))                            :: z,z2
  REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: dr_vct
  REAL(KIND(0D0))                            :: dr,dpi
  REAL(KIND(0E0)),DIMENSION(:),ALLOCATABLE   :: r_vct
  REAL(KIND(0E0))                            :: r,pi
  INTEGER                              :: i,n


  n = 100
  ALLOCATE (dz_vct(n),z_vct(n),dr_vct(n),r_vct(n))
! REAL 
  print*,'REAL arithmetics'
  print*,'KIND(r):',KIND(r)
  print*,'DIGITS(r):',DIGITS(r)
  print*,'EPSILON(r):',EPSILON(r)
  print*
  pi = ACOS(-1.)
  print*,'pi = ACOS(-1.):',pi
  pi = ACOS(-1d0)
  print*,'pi = ACOS(-1d0):',pi
  r = SQRT(pi)
  print*,'r = SQRT(pi):',r
  print*
  r_vct = pi
  print*,'KIND(r_vct):',KIND(r_vct)
  print*,'SHAPE(r_vct):',SHAPE(r_vct)
  print*,'r_vct(:)=',pi
  r=0d0
  DO i=1,n
     r = r + r_vct(i)*r_vct(i)
  END DO
  print*,'Inner product of r_vct via do loop',r
  print*,'and via DOT_PRODUCT(r_vct,r_vct)',DOT_PRODUCT(r_vct,r_vct)
  print*
! DOUBLE PRECISION
  print*,'DOUBLE PRECISION arithmetics'
  print*,'KIND(dr):',KIND(dr)
  print*,'EPSILON(dr):',EPSILON(dr)
  print*
  dpi = ACOS(-1d0)
  print*,'dpi = ACOS(-1d0):',dpi
  dpi = ACOS(-1.)
  print*,'dpi = ACOS(-1.):',dpi
  dpi = DACOS(-1d0)
  print*,'dpi = DACOS(-1d0):',dpi
  dr = SQRT(dpi)
  print*,'dr = SQRT(dpi):',dr
  dr = DSQRT(dpi)
  print*,'dr = DSQRT(dpi):',dr
  print*
  dr_vct = dpi
  print*,'KIND(dr_vct):',KIND(dr_vct)
  print*,'SHAPE(dr_vct):',SHAPE(dr_vct)
  print*,'dr_vct(:)=',dpi
  dr=0d0
  DO i=1,n
     dr = dr + dr_vct(i)*dr_vct(i)
  END DO
  print*,'Inner product of dr_vct via do loop',dr
  print*,'and via DOT_PRODUCT(dr_vct,dr_vct)',DOT_PRODUCT(dr_vct,dr_vct)
  print*
! COMPLEX
  print*,'COMPLEX arithmetics'
  print*,'KIND(z):',KIND(z)
  print*
  z = CMPLX(1.0,pi)
  print*,'z = CMPLX(1.0,pi):',z
  z = DCMPLX(1.0,pi)
  print*,'z = DCMPLX(1.0,pi):',z
  z = CMPLX(pi,pi)
  print*,'z = CMPLX(pi,pi):',z
  z = DCMPLX(pi,pi)
  print*,'z = DCMPLX(pi,pi):',z
  z2 = SQRT(z)
  print*,'z2 = SQRT(z):',z2
  z2 = EXP(z)
  print*,'z2 = EXP(z):',z2
  z2 = EXP(DCMPLX(pi,pi))
  print*,'z2 = EXP(DCMPLX(pi,pi)):',z2
  z2 = CDEXP(DCMPLX(pi,pi))
  print*,'z2 = CDEXP(DCMPLX(pi,pi)):',z2
  z2 = CONJG(z)*z
  print*,'z2 = CONJG(z)*z:',z2
  print*
  print*,'KIND(z_vct):',KIND(z_vct)
  print*,'SHAPE(z_vct):',SHAPE(z_vct)
  z_vct = CMPLX(pi,1.0)
  print*,'z_vct(:)=',z_vct(1)
  r=0.0
  DO i=1,n
     r = r + REAL(CONJG(z_vct(i))*z_vct(i))
  END DO
  print*,'Inner product of dz_vct via do loop',r
  print*,'and via DOT_PRODUCT(z_vct,z_vct)',&
       DOT_PRODUCT(z_vct,z_vct)
  print*,'and via DOT_PRODUCT(CONJG(dz_vct),dz_vct)',&
       DOT_PRODUCT(CONJG(z_vct),z_vct)
  print*,'ABS(DOT_PRODUCT(CONJG(z_vct),z_vct))',&
       ABS(DOT_PRODUCT(CONJG(z_vct),z_vct))
  print*
  print*,'DOUBLE PRECISION COMPLEX arithmetics'
  print*,'KIND(dz):',KIND(dz)
!  print*,'DIGITS(dz):',DIGITS(dz)
  print*
  dz = CMPLX(1.0,pi)
  print*,'dz = CMPLX(1.0,pi):',dz
  dz = DCMPLX(1.0,pi)
  print*,'dz = DCMPLX(1.0,pi):',dz
  dz = CMPLX(pi,pi)
  print*,'dz = CMPLX(pi,pi):',dz
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
  dz2 = CONJG(dz)*dz
  print*,'dz = CONJG(dz)*dz:',dz2
  print*
  print*,'KIND(dz_vct):',KIND(dz_vct)
  print*,'SHAPE(dz_vct):',SHAPE(dz_vct)
  dz_vct = DCMPLX(pi,1d0)
  print*,'dz_vct(:)=',dz_vct(1)
  dr=0d0
  DO i=1,n
     dr = dr + DBLE(CONJG(dz_vct(i))*dz_vct(i))
  END DO
  print*,'Inner product of dz_vct via do loop',dr
  print*,'and via DOT_PRODUCT(dz_vct,dz_vct)',&
       DOT_PRODUCT(dz_vct,dz_vct)
  print*,'ABS(DOT_PRODUCT(dz_vct,dz_vct))',&
       ABS(DOT_PRODUCT(dz_vct,dz_vct))
  print*,'and via DOT_PRODUCT(CONJG(dz_vct),dz_vct)',&
       DOT_PRODUCT(CONJG(dz_vct),dz_vct)
  print*,'ABS(DOT_PRODUCT(CONJG(dz_vct),dz_vct))',&
       ABS(DOT_PRODUCT(CONJG(dz_vct),dz_vct))
  
END PROGRAM main

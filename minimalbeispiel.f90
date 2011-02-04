PROGRAM main

  USE tic_toc! counts calculation time

  IMPLICIT none

  COMPLEX(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: dz_vct
  COMPLEX(KIND(0D0))                            :: dz,dz2
  COMPLEX(KIND(0E0)),DIMENSION(:),ALLOCATABLE   :: z_vct
  COMPLEX(KIND(0E0))                            :: z,z2
  REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: dr_vct
  REAL(KIND(0D0))                            :: dr,dpi
  REAL(KIND(0E0)),DIMENSION(:),ALLOCATABLE   :: r_vct,r_vct2
  REAL(KIND(0E0))                            :: r,pi
  INTEGER                              :: i,n,m,j
  INTEGER                              :: c1
  CHARACTER (256)                       :: csz

  n = 10 ! vector size
  m = 1 ! loop counts
  ALLOCATE (dz_vct(n),z_vct(n),dr_vct(n),r_vct(n),r_vct2(n))
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
  print*,'r_vct(:)=',r_vct(1)
  csz = 'DO LOOP::'
  CALL TIC(c1)
  DO j=1,m
     r=0d0
     DO i=1,n
        r = r + r_vct(i)*r_vct(i)
     END DO
  END DO
  print*,'Inner product of r_vct via do loop',r
  CALL TOC(c1,csz)
  DO i=1,n
     r_vct2(i) = REAL(i)
  END DO
  print*,'r_vct2:',r_vct2
  r_vct = r_vct * r_vct2
  print*,'r_vct = r_vct * r_vct',r_vct(:)
  STOP
  csz = 'DOT_PROD::'
  CALL TIC(c1)
  DO j=1,m
     r = DOT_PRODUCT(r_vct,r_vct)
  END DO
  print*,'and via DOT_PRODUCT(r_vct,r_vct)',r
  CALL TOC(c1,csz)
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
  dr_vct = 1d0/dpi
  print*,'KIND(dr_vct):',KIND(dr_vct)
  print*,'SHAPE(dr_vct):',SHAPE(dr_vct)
  print*,'dr_vct(:)=',dr_vct(1)
  csz = 'DO LOOP::'
  CALL TIC(c1)
  DO j=1,m
     dr=0d0
     DO i=1,n
        dr = dr + dr_vct(i)*dr_vct(i)
     END DO
  END DO
  print*,'Inner product of dr_vct via do loop',dr
  CALL TOC(c1,csz)
  csz = 'DOT_PROD::'
  CALL TIC(c1)
  DO j=1,m
     dr = DOT_PRODUCT(dr_vct,dr_vct)
  END DO
  print*,'and via DOT_PRODUCT(dr_vct,dr_vct)',dr
  CALL TOC(c1,csz)
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
  z_vct = 1/CMPLX(pi,1.0)
  print*,'z_vct(:)=',z_vct(1)
  csz = 'DO LOOP::'
  CALL TIC(c1)
  DO j=1,m
     r=0.
     DO i=1,n
        r = r + REAL(CONJG(z_vct(i))*z_vct(i))
     END DO
  END DO
  print*,'Inner product of dz_vct via do loop',r
  CALL TOC(c1,csz)
  csz = 'DOT_PROD::'
  CALL TIC(c1)
  DO j=1,m
     r = DOT_PRODUCT(z_vct,z_vct)
  END DO
  print*,'and via DOT_PRODUCT(z_vct,z_vct)',r
  CALL TOC(c1,csz)

  print*
  print*,'DOUBLE PRECISION COMPLEX arithmetics'
  print*,'KIND(dz):',KIND(dz)
  print*,'KIND(pi):',KIND(pi)
  print*,'KIND(dpi):',KIND(dpi)
  print*,'pi:',pi
  print*,'dpi:',dpi

  !  print*,'DIGITS(dz):',DIGITS(dz)
  print*
  dz = CMPLX(1.0,dpi)
  print*,'dz = CMPLX(1.0,dpi):',dz
  dz = CMPLX(1.0,pi)
  print*,'dz = DCMPLX(1.0,pi):',dz
  dz = CMPLX(dpi,pi)
  print*,'dz = CMPLX(dpi,pi):',dz
  dz = DCMPLX(dpi,pi)
  print*,'dz = DCMPLX(dpi,pi):',dz
  dz = CMPLX(dpi,dpi)
  print*,'dz = CMPLX(dpi,dpi):',dz
  dz = DCMPLX(dpi,dpi)
  print*,'dz = DCMPLX(dpi,dpi):',dz
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
  dz2 = EXP(DCMPLX(dpi,dpi))
  print*,'dz = EXP(DCMPLX(dpi,dpi)):',dz2
  dz2 = CDEXP(DCMPLX(dpi,dpi))
  print*,'dz = CDEXP(DCMPLX(dpi,dpi)):',dz2
  dz2 = CONJG(dz)*dz
  print*,'dz = CONJG(dz)*dz:',dz2
  print*
  print*,'KIND(dz_vct):',KIND(dz_vct)
  print*,'SHAPE(dz_vct):',SHAPE(dz_vct)
  dz_vct = 1D0/DCMPLX(pi,1d0)
  print*,'dz_vct(:)=',dz_vct(1)

  csz = 'DO LOOP::'
  CALL TIC(c1)
  DO j=1,m
     dr=0d0
     DO i=1,n
        dr = dr + DBLE(CONJG(dz_vct(i))*dz_vct(i))
     END DO
  END DO
  print*,'Inner product of dz_vct via do loop',dr
  CALL TOC(c1,csz)
  csz = 'DOT_PROD::'
  CALL TIC(c1)
  DO j=1,m
     dr = DOT_PRODUCT(dz_vct,dz_vct)
  END DO
  print*,'and via DOT_PRODUCT(dz_vct,dz_vct)',dr
  CALL TOC(c1,csz)

  print*,'COS(ACOS(1./3.)),COS(ACOS(1d0/3))',COS(ACOS(1./3.)),&
       COS(ACOS(1d0/3))

END PROGRAM main

!!!$ Gauss-Legendre "integration"
!!!$ Numerical Recipies (1986), p. 125, subroute GAULEG
!!!$ Returns the abcissas (x) and weights(w) of length n
!!!$ For a gauss-quadrature integration.
!!!$ Provideds x_i and w_i for the formula: 
!!!$ \sum_{x_1}^{x_2} f(x) dx = \sum_{i=1}^N w_i f(x_i)
SUBROUTINE gauleg(x1,x2,x,w,n)
use alloci, only: prec
  IMPLICIT none
  INTEGER (KIND = 4)  :: n
  REAL (prec)    :: x1,x2,x(n),w(n)
  REAL (prec)    :: EPS
  PARAMETER (EPS=3.d-14)
  INTEGER (KIND = 4)  :: i,j,m
  REAL (prec)    :: p1,p2,p3,pp,xl,xm,z,z1
  m=(n+1)/2
  xm=0.5d0*(x2+x1)
  xl=0.5d0*(x2-x1)
  do i=1,m
     z=cos(3.141592653589793d0*(REAL(i)-.25d0)/(REAL(n)+.5d0))
1    continue
     p1=1.d0
     p2=0.d0
     do j=1,n
        p3=p2
        p2=p1
        p1=((2.d0*REAL(j)-1.d0)*z*p2-(REAL(j)-1.d0)*p3)/REAL(j)
     END DO
     pp=REAL(n)*(z*p1-p2)/(z*z-1.d0)
     z1=z
     z=z1-p1/pp
     if(ABS(z-z1).gt.EPS)goto 1
     x(i)=xm-xl*z
     x(n+1-i)=xm+xl*z
     w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
     w(n+1-i)=w(i)
  END DO
  return
END SUBROUTINE gauleg

!!!$ Gauss-Laguerre "Integration": Computation of abcissas
!!!$ and weights for the numerical integration of the k-integral
!!!$ for the back-transformation to the y-space. From
!!!$ Numerical Recipies, 1992, p. 146
SUBROUTINE gaulag(x,w,n,alf)
use alloci, only: prec
  IMPLICIT none
  INTEGER (KIND = 4) ::n,MAXIT
  REAL (prec) :: alf,w(n),x(n)
  REAL (prec) :: EPS
  PARAMETER (EPS=3.e-14,MAXIT=10)
  INTEGER (KIND = 4) :: i,its,j
  REAL (prec) :: ai,gammln
  REAL (prec) :: p1,p2,p3,pp,z,z1
  do i=1,n
     if(i.eq.1)then
        z=(1.+alf)*(3.+.92*alf)/(1.+2.4*n+1.8*alf)
     else if(i.eq.2)then
        z=z+(15.+6.25*alf)/(1.+.9*alf+2.5*n)
     else
        ai=REAL(i-2)
        z=z+((1.+2.55*ai)/(1.9*ai)+1.26*ai*alf/ &
             (1.+3.5*ai))*(z-x(i-2))/(1.+.3*alf)
     endif
     do its=1,MAXIT
        p1=1.
        p2=0.
        do j=1,n
           p3=p2
           p2=p1
           p1=((REAL(2*j-1)+alf-z)*p2-(REAL(j-1)+alf)*p3)/REAL(j)
        END do
        pp=(REAL(n)*p1-(REAL(n)+alf)*p2)/z
        z1=z
        z=z1-p1/pp
        if(ABS(z-z1).le.EPS)goto 1
     END do
     print*,'too many iterations in gaulag'
1    x(i)=z
     w(i)=-EXP(gammln(alf+REAL(n,prec))-gammln(real(n,prec)))/(pp*REAL(n)*p2)
  END do
  return
END SUBROUTINE gaulag

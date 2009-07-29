      SUBROUTINE gauleg(x1,x2,x,w,n)
      integer * 4 n
      real    * 8 x1,x2,x(n),w(n)
      real    * 8 EPS
      PARAMETER (EPS=3.d-14)
      integer * 4 i,j,m
      real    * 8 p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592653589793d0*(dble(i)-.25d0)/(dble(n)+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*dble(j)-1.d0)*z*p2-(dble(j)-1.d0)*p3)/dble(j)
11        continue
          pp=dble(n)*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(dabs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END

FUNCTION gammln(xx)
use alloci, only: prec
  REAL (prec)  :: xx,gammln
  INTEGER (KIND=4)  :: j
  REAL (prec) :: ser,stp,tmp,x,y,cof(6)
  SAVE cof,stp
  DATA cof,stp/76.18009172947146,-86.50532032941677, &
       24.01409824083091,-1.231739572450155,.1208650973866179e-2, &
       -.5395239384953e-5,2.5066282746310005/
  x=xx
  y=x
  tmp=x+5.5_prec
  tmp=(x+0.5)*LOG(tmp)-tmp
  ser=1.000000000190015_prec
  do j=1,6
     y=y+1.
     ser=ser+cof(j)/y
  END DO
  gammln=tmp+LOG(stp*ser/x)
  return
END FUNCTION gammln

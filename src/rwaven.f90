subroutine rwaven()
! compute wave numbers 
!
! Andreas Kemna                                            20-Dec-1993
! last change                                              19-Jun-1998
! ..............................................................................
  use alloci, only: prec
  USE datmod
  USE electrmod
  USE elemmod
  USE wavenmod
  USE errmod

  IMPLICIT none

! ..............................................................................
! internals
! index vars
  INTEGER           i, j
! electrode numbers
  INTEGER           elec( 4 )
! dummy vars
  INTEGER           ganz, lanz
  REAL (prec)       kwn0, dum, xke( 4 ), yke( 4 ), x21, y21

! ..............................................................................
  amin = 1d9
  amax = 0d0
  do i=1,nanz
! get current and potential electrodes
     elec(1) = mod(strnr(i),10000)
     elec(2) = (strnr(i)-elec(1))/10000
     elec(3) = mod(vnr(i),10000)
     elec(4) = (vnr(i)-elec(3))/10000
! get electrode distances
     do j=1,4
        if (elec(j).gt.0) xke(j) = sx(snr(enr(elec(j))))
        if (elec(j).gt.0) yke(j) = sy(snr(enr(elec(j))))
     end do
     if (elec(3).gt.0) then
        do j=1,2
           if (elec(j).gt.0) then
              x21  = xke(3)-xke(j)
              y21  = yke(3)-yke(j)
              dum  = SQRT(x21*x21+y21*y21)
              amin = INT(MIN(dum,amin)) ! why fuckin INT here?
              amax = INT(MAX(dum,amax))
           end if
        end do
     end if
     if (elec(4).gt.0) then
        do j=1,2
           if (elec(j).gt.0) then
              x21  = xke(4)-xke(j)
              y21  = yke(4)-yke(j)
              dum  = SQRT(x21*x21+y21*y21)
              amin = INT(MIN(dum,amin))
              amax = INT(MAX(dum,amax))
           end if
        end do
     end if
  end do

! errors
  if (amin.lt.1d-12.or.amax.lt.1d-12) then
     WRITE (fetxt,*)'check for duplicate electrodes'
     errnr = 73
     goto 1000
  end if

! compute wave numbers
! Andreas Kemna Dissertation pp.163
  lanz   = 4
! number of abcissas for Gauss Legende integration (AK Diss pp.163)
  ganz   = int(real(6d0*LOG10(amax/amin)))
! need at least 2 abcissas
  ganz   = max0(ganz,2)
  ganz   = min0(ganz,30 - lanz)
! total k number
  kwnanz = ganz+lanz
! k_0, AK Diss pp. 163
  kwn0   = 1d0/(2d0*amin)
  ALLOCATE (kwn(kwnanz),kwnwi(kwnanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation kwn'
     errnr = 94
     RETURN
  END IF

! Gauss integeration
! variable substitution in order to overcome the singularity at zero
! k' = (k / k_0)^{1/2}
! integration range changes from (0, k_0) to (0, 1)
  call gauleg(0._prec,1._prec,kwn(1),kwnwi(1),ganz)

! see AK Diss, pp. 161
! weighting factor also needs substitution: w_n = 2 k_0 k'_n w'_n / pi
! the pi vanishes when evaluating the whole integral -> no need to use it
! anywhere.
! k_n = k_0 k'_n^2
  do i=1,ganz
     kwnwi(i) = 2._prec*kwn0*kwnwi(i)*kwn(i)
     kwn(i)   = kwn0*kwn(i)*kwn(i)
  end do

! Laguerre integration
! compute abcissas and weights for the numerical integration
! of the upper integral (see AK Diss, p. 161ff.)
! the alpha = 0d0 value removes the exp function in the
! Gauss-Laguere formula (Numerical recipies, p. 144, 1992)
  call gaulag(kwn((ganz+1):kwnanz),kwnwi((ganz+1):kwnanz),lanz,0._prec)

  do i=ganz+1,ganz+lanz
     kwnwi(i) = kwn0*kwnwi(i)*EXP(kwn(i))
     kwn(i)   = kwn0*(kwn(i)+1._prec)
  end do
  errnr = 0
  return
! ..............................................................................
! error messages

1000 return

end subroutine rwaven

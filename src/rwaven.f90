subroutine rwaven()

!!!$     Unterprogramm zur Berechnung der Wellenzahlwerte.

!!!$     Andreas Kemna                                            20-Dec-1993
!!!$     Letzte Aenderung   19-Jun-1998

!!!$.....................................................................

  USE datmod
  USE electrmod
  USE elemmod
  USE wavenmod
  USE errmod

  IMPLICIT none


!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Indexvariablen
  INTEGER (KIND=4) ::  i,j

!!!$     Elektrodennummern
  INTEGER (KIND=4) ::     elec(4)

!!!$     Hilfsvariablen
  INTEGER (KIND=4) ::     ganz,lanz
  REAL (KIND(0D0)) ::    kwn0,dum,xke(4),yke(4),x21,y21

!!!$.....................................................................

!!!$     'amin', 'amax' bestimmen
  amin = 1d9
  amax = 0d0

  do i=1,nanz

!!!$     Stromelektroden bestimmen
     elec(1) = mod(strnr(i),10000)
     elec(2) = (strnr(i)-elec(1))/10000

!!!$     Messelektroden bestimmen
     elec(3) = mod(vnr(i),10000)
     elec(4) = (vnr(i)-elec(3))/10000

!!!$     Abstaende bestimmen
     do j=1,4
        if (elec(j).gt.0) xke(j) = sx(snr(enr(elec(j))))
        if (elec(j).gt.0) yke(j) = sy(snr(enr(elec(j))))
     end do

     if (elec(3).gt.0) then
        do j=1,2
           if (elec(j).gt.0) then
              x21  = xke(3)-xke(j)
              y21  = yke(3)-yke(j)
              dum  = dsqrt(x21*x21+y21*y21)
              amin = dmin1(dum,amin)
              amax = dmax1(dum,amax)
           end if
        end do
     end if

     if (elec(4).gt.0) then
        do j=1,2
           if (elec(j).gt.0) then
              x21  = xke(4)-xke(j)
              y21  = yke(4)-yke(j)
              dum  = dsqrt(x21*x21+y21*y21)
              amin = dmin1(dum,amin)
              amax = dmax1(dum,amax)
           end if
        end do
     end if
  end do

!!!$     Ggf. Fehlermeldung
  if (amin.lt.1d-12.or.amax.lt.1d-12) then
     fetxt = ' '
     errnr = 73
     goto 1000
  end if

!!!$     Wellenzahlen bestimmen
!!!$ AK Diss p. 163
  lanz   = 4
!!!$     ak        ganz   = int(real(6d0*dlog10(amax)))
!!!$ Number of abcissas for Gauss-Legende-Integration (AK Diss p. 163)
  ganz   = int(real(6d0*dlog10(amax/amin)))
!!!$ Need at least 2 abcissas
  ganz   = max0(ganz,2)
  ganz   = min0(ganz,30 - lanz)
!!!$ Overall k number
  kwnanz = ganz+lanz
!!!$ k_0, AK Diss p. 163
  kwn0   = 1d0/(2d0*amin)

  ALLOCATE (kwn(kwnanz),kwnwi(kwnanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation kwn'
     errnr = 94
     RETURN
  END IF

!!!$ REMEMBER: We only compute the k values and the weighting factors
!!!$ here. No actual integration is done here as we still need to caluclate 
!!!$ the potentials for the k-values.

!!!$     Gauss-Integeration
!!!$ Variable substitution in order to overcome the singularity at zero
!!!$ k' = (k / k_0)^{1/2}
!!!$ Integration range changes from (0, k_0) to (0, 1)
  call gauleg(0d0,1d0,kwn(1),kwnwi(1),ganz)

!!!$ See Ak-Diss(2000), p. 161
!!!$ weighting factor also needs substitution: w_n = 2 k_0 k'_n w'_n / pi
!!!$ The pi vanishes when evaluating the whole integral -> no need to use it
!!!$ anywhere.
!!!$ k_n = k_0 k'_n^2
  do i=1,ganz
     kwnwi(i) = 2d0*kwn0*kwnwi(i)*kwn(i)
     kwn(i)   = kwn0*kwn(i)*kwn(i)
  end do

!!!$     Laguerre-Integration
!!!$ Compute abcissas and weights for the numerical integration
!!!$ of the upper integral (see AK Diss, p. 161ff.)
!!!$ the alpha = 0d0 value removes the exp function in the
!!!$ Gauss-Laguere formula (Numerical recipies, p. 144, 1992)
  call gaulag(kwn(ganz+1),kwnwi(ganz+1),lanz,0d0)

  do i=ganz+1,ganz+lanz
     kwnwi(i) = kwn0*kwnwi(i)*dexp(kwn(i))
     kwn(i)   = kwn0*(kwn(i)+1d0)
  end do

  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

1000 return

end subroutine rwaven

      subroutine rwaven()
      
c     Unterprogramm zur Berechnung der Wellenzahlwerte.

c     Andreas Kemna                                            20-Dec-1993
c     Letzte Aenderung   19-Jun-1998

c.....................................................................

      USE datmod
      USE electrmod
      USE elemmod
      USE wavenmod
      USE errmod

      IMPLICIT none


c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Indexvariablen
      integer         * 4     i,j

c     Elektrodennummern
      integer         * 4     elec(4)

c     Hilfsvariablen
      integer		  * 4     ganz,lanz
      real            * 8     kwn0,dum,
     1     xke(4),yke(4),
     1     x21,y21

c.....................................................................

c     'amin', 'amax' bestimmen
      amin = 1d9
      amax = 0d0

      do i=1,nanz

c     Stromelektroden bestimmen
         elec(1) = mod(strnr(i),10000)
         elec(2) = (strnr(i)-elec(1))/10000

c     Messelektroden bestimmen
         elec(3) = mod(vnr(i),10000)
         elec(4) = (vnr(i)-elec(3))/10000

c     Abstaende bestimmen
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

c     Ggf. Fehlermeldung
      if (amin.lt.1d-12.or.amax.lt.1d-12) then
         fetxt = ' '
         errnr = 73
         goto 1000
      end if

c     Wellenzahlen bestimmen
      lanz   = 4
c     ak        ganz   = int(real(6d0*dlog10(amax)))
      ganz   = int(real(6d0*dlog10(amax/amin)))
      ganz   = max0(ganz,2)
      kwnanz = ganz+lanz
      kwn0   = 1d0/(2d0*amin)

      ALLOCATE (kwn(kwnanz),kwnwi(kwnanz),stat=errnr)
      IF (errnr /= 0) THEN
         fetxt = 'Error memory allocation kwn'
         errnr = 94
         RETURN
      END IF
 
c     Gauss-Integeration
      call gauleg(0d0,1d0,kwn(1),kwnwi(1),ganz)
      
      do i=1,ganz
         kwnwi(i) = 2d0*kwn0*kwnwi(i)*kwn(i)
         kwn(i)   = kwn0*kwn(i)*kwn(i)
      end do

c     Laguerre-Integration
      call gaulag(kwn(ganz+1),kwnwi(ganz+1),lanz,0d0)

      do i=ganz+1,ganz+lanz
         kwnwi(i) = kwn0*kwnwi(i)*dexp(kwn(i))
         kwn(i)   = kwn0*(kwn(i)+1d0)
      end do

      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 1000 return

      end

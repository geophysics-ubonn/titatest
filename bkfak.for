      subroutine bkfak()

c     Unterprogramm zur Berechnung der Konfigurationsfaktoren
c     (bei Standard-Geometrie!).

c     Andreas Kemna                                            02-May-1996
c     Letzte Aenderung   29-Apr-2003

c.....................................................................
      USE datmod

      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'electr.fin'

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Indexvariable
      integer         * 4     i,fp

c     Elektrodennummern
      integer         * 4     elec1,elec2,
     1     elec3,elec4

c     Koordinaten
      real            * 8     xk(4),yk(4)

c     Pi
      real            * 8     pi

c     Hilfsvariablen
      real            * 8     dx,dym,dyp,dum,
     1     r1,r2,r3,r4

c.....................................................................

      pi = dacos(-1d0)

      CALL get_unit(fp)
      OPEN (fp,FILE='tmp.kfak',STATUS='replace')
      
      do i=1,nanz
         elec1 = mod(strnr(i),10000)
         elec2 = (strnr(i)-elec1)/10000

         elec3 = mod(vnr(i),10000)
         elec4 = (vnr(i)-elec3)/10000

         if (elec1.gt.0) then
            xk(1) = sx(snr(enr(elec1)))
            yk(1) = sy(snr(enr(elec1)))
         end if

         if (elec2.gt.0) then
            xk(2) = sx(snr(enr(elec2)))
            yk(2) = sy(snr(enr(elec2)))
         end if

         if (elec3.gt.0) then
            xk(3) = sx(snr(enr(elec3)))
            yk(3) = sy(snr(enr(elec3)))
         end if

         if (elec4.gt.0) then
            xk(4) = sx(snr(enr(elec4)))
            yk(4) = sy(snr(enr(elec4)))
         end if

         if (elec3.gt.0.and.elec1.gt.0) then
            dx  = xk(3)-xk(1)
            dym = yk(3)-yk(1)
            dyp = yk(3)+yk(1)
            dym = 1d0/dsqrt(dx*dx+dym*dym)
            dyp = 1d0/dsqrt(dx*dx+dyp*dyp)
            r4  = dym+dyp
         else
            r4  = 0d0
         end if

         if (elec3.gt.0.and.elec2.gt.0) then
            dx  = xk(3)-xk(2)
            dym = yk(3)-yk(2)
            dyp = yk(3)+yk(2)
            dym = 1d0/dsqrt(dx*dx+dym*dym)
            dyp = 1d0/dsqrt(dx*dx+dyp*dyp)
            r3  = dym+dyp
         else
            r3  = 0d0
         end if

         if (elec4.gt.0.and.elec1.gt.0) then
            dx  = xk(4)-xk(1)
            dym = yk(4)-yk(1)
            dyp = yk(4)+yk(1)
            dym = 1d0/dsqrt(dx*dx+dym*dym)
            dyp = 1d0/dsqrt(dx*dx+dyp*dyp)
            r2  = dym+dyp
         else
            r2  = 0d0
         end if

         if (elec4.gt.0.and.elec2.gt.0) then
            dx  = xk(4)-xk(2)
            dym = yk(4)-yk(2)
            dyp = yk(4)+yk(2)
            dym = 1d0/dsqrt(dx*dx+dym*dym)
            dyp = 1d0/dsqrt(dx*dx+dyp*dyp)
            r1  = dym+dyp
         else
            r1  = 0d0
         end if

         dum = (r1-r2) - (r3-r4)

         if (dabs(dum).eq.0d0) then
            fetxt = 'index '
            write(fetxt(7:10),'(i4)') i
            errnr = 93
            goto 1000
c     tmp                write(*,'(i4)') i
c     tmp                dum = 1d-12
         end if

         kfak(i) = 4d0*pi / dum
         
      end do
      WRITE (fp,'(I8,2X,G12.4)')(i,kfak(i),i=1,nanz)
      CLOSE (fp)

      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 1000 return

      end

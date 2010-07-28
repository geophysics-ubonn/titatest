      subroutine bsmatm()

c     Unterprogramm belegt die Rauhigkeitsmatrix.

c     Andreas Kemna                                            29-Feb-1996
c     Letzte Aenderung   04-Nov-2008

c.....................................................................

      USE alloci
      USE datmod
      USE invmod
      USE modelmod

      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'konv.fin'

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Variablen zur Beruecksichtigung von Diskontinuitaeten (keine
c     Glaettung in x bzw. z-Richtung)
      integer         * 4     ndis_z,idis_z(3),m,
     1     ndis_x,idis_x(4)
      logical         * 4     lup,ldown,lleft,lright
      real            * 8     alfdis

c     Hilfsvariablen
      real            * 8     dum,
     1     dzleft,dzright,
     1     xleft,xmean,xright,
     1     dxup,dxdown,
     1     zup,zmean,zdown

      integer         * 4     i,j,l
      
c     Hilfsfunction
      integer         * 4     k

      k(i,j) = (i-1) * nx + j

c.....................................................................

      IF (.NOT.ALLOCATED (smatm)) ALLOCATE (smatm(manz,3))
      ndis_z = 0
c     ak BAW
c     ak        ndis_z    = 2
      idis_z(1) = 15
      idis_z(2) = 18

      ndis_x = 0
c     ak Bohrloch-Effekt
c     ak        ndis_x    = 4
      idis_x(1) = 4
      idis_x(2) = 6
      idis_x(3) = nx-4
      idis_x(4) = nx-2

c     ak fuer Christoph (Wald)
c     ak        ndis_z    = 3
      idis_z(1) = 6
      idis_z(2) = 12
      idis_z(3) = 19
      
      IF (.NOT.ALLOCATED (smatm)) ALLOCATE (smatm(manz,3))

c     Rauhigkeitsmatrix auf Null setzen
      do i=1,manz
         do j=1,3
            smatm(i,j) = 0d0
         end do
      end do

      do i=1,nz
         lup   = .true.
         ldown = .true.

         do m=1,ndis_z
            if (i  .eq.idis_z(m)) lup  =.false.
            if (i+1.eq.idis_z(m)) ldown=.false.
         end do

         do j=1,nx
            lleft  = .true.
            lright = .true.

            do m=1,ndis_x
               if (j  .eq.idis_x(m)) lleft =.false.
               if (j+1.eq.idis_x(m)) lright=.false.
            end do

c     Beitrag von Wx^t*Wx zur Rauhigkeitsmatrix
            dzleft  = dabs( sy(snr(nrel(k(i,j),4)))
     1           -sy(snr(nrel(k(i,j),1))))
            dzright = dabs( sy(snr(nrel(k(i,j),3)))
     1           -sy(snr(nrel(k(i,j),2))))

            xmean = 0d0
            do l=1,4
               xmean = xmean + sx(snr(nrel(k(i,j),l)))
            end do
            xmean = xmean/4d0


            if (j.gt.1) then
               if (lleft) then
                  alfdis = 1d0
               else
c     ak
                  alfdis = 1d-3
               end if

               xleft = 0d0
               do l=1,4
                  xleft = xleft + sx(snr(nrel(k(i,j-1),l)))
               end do

               xleft = xleft/4d0

               smatm(k(i,j),1) = alfdis*alfx *
     1              dzleft/dabs(xmean-xleft)
            end if

            if (j.lt.nx) then
               if (lright) then
                  alfdis = 1d0
               else
c     ak
                  alfdis = 1d-3
               end if

               xright = 0d0
               do l=1,4
                  xright = xright + sx(snr(nrel(k(i,j+1),l)))
               end do

               xright = xright/4d0
               dum    = alfdis*alfx * dzright/dabs(xright-xmean)

               smatm(k(i,j),1) = smatm(k(i,j),1) + dum
               smatm(k(i,j),2) = -dum
            end if

c     Beitrag von Wz^t*Wz zur Rauhigkeitsmatrix
            dxup   = dabs( sx(snr(nrel(k(i,j),3)))
     1           -sx(snr(nrel(k(i,j),4))))
            dxdown = dabs( sx(snr(nrel(k(i,j),2)))
     1           -sx(snr(nrel(k(i,j),1))))

            zmean = 0d0
            do l=1,4
               zmean = zmean + sy(snr(nrel(k(i,j),l)))
            end do
            zmean = zmean/4d0

            if (i.gt.1) then
               if (lup) then
                  alfdis = 1d0
               else
                  alfdis = 0d0
               end if

               zup = 0d0
               do l=1,4
                  zup = zup + sy(snr(nrel(k(i-1,j),l)))
               end do

               zup = zup/4d0
               dum = alfdis*alfz * dxup/dabs(zup-zmean)

               smatm(k(i,j),1) = smatm(k(i,j),1) + dum
            end if

            if (i.lt.nz) then
               if (ldown) then
                  alfdis = 1d0
               else
                  alfdis = 0d0
               end if

               zdown = 0d0
               do l=1,4
                  zdown = zdown + sy(snr(nrel(k(i+1,j),l)))
               end do

               zdown = zdown/4d0
               dum   = alfdis*alfz * dxdown/dabs(zmean-zdown)

               smatm(k(i,j),1) = smatm(k(i,j),1) + dum
               smatm(k(i,j),3) = -dum
            end if

         end do
      end do

      return
      end

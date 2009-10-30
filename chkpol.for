      subroutine chkpol(lsetup)

c     Unterprogramm zum automatischen Korrigieren der Polaritaeten.

c     Andreas Kemna                                            07-Dec-1996
c     Letzte Aenderung   07-Nov-1997
      
c.....................................................................
      IMPLICIT none
      INCLUDE 'err.fin'
      INCLUDE 'path.fin'
      INCLUDE 'parmax.fin'
      INCLUDE 'dat.fin'
      INCLUDE 'inv.fin'

c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Hilfsschalter
      logical         * 4     lsetup

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Messelektrodennummern
      integer         * 4     elec3,elec4

c     Hilfsvariablen
      integer         * 4     i,idat,isig

c     Real-, Imaginaerteile
      real            * 8     redat,imdat,
     1     resig,imsig

c     Pi
      real            * 8     pi

c.....................................................................

      pi = dacos(-1d0)

c     'run.ctr' oeffnen
      open(10,file=ramd(1:lnramd)//slash(1:1)//'run.ctr',status='old')
 1    read(10,*,end=2)
      goto 1
 2    backspace(10)

      do i=1,nanz
         wdfak(i) = 1

c     Messelektroden bestimmen
         elec3 = mod(vnr(i),10000)
         elec4 = (vnr(i)-elec3)/10000

c     Logarithmierte Betraege in den Realteilen,
c     Phasen (in rad) in den Imaginaerteilen
         redat = dble(dat(i))
         imdat = dimag(dat(i))
         resig = dble(sigmaa(i))
         imsig = dimag(sigmaa(i))

c     Phasenbereich checken
         if (imdat.gt.pi/2d0) then
            idat = -1
         else if (imdat.le.-pi/2d0) then
            idat = 1
         else
            idat = 0
         end if

         if (imsig.gt.pi/2d0) then
            isig = -1
         else if (imsig.le.-pi/2d0) then
            isig = 1
         else
            isig = 0
         end if

         if (idat.eq.0.and.isig.ne.0) then

c     Falls lpol=.true., angenommene Polaritaet des Messdatums falsch,
c     ggf. Korrektur; auf jeden Fall Polaritaetswechsel
            vnr(i)    = elec3*10000 + elec4
            imsig     = imsig + dble(isig)*pi
            sigmaa(i) = dcmplx(resig,imsig)
            volt(i)   = cdexp(-sigmaa(i))
            
            if (lpol) then
               write(10,'(i4,a30)')
     1              i,' : correct and change polarity'
               if (.not.lsetup) wdfak(i)=0
            else
               imdat  = imdat - dsign(pi,imdat)
               dat(i) = dcmplx(redat,imdat)

               write(10,'(i4,a18)')
     1              i,' : change polarity'
               wdfak(i) = 0
            end if

         else if (idat.ne.0.and.isig.eq.0) then

c     Falls lpol=.true., angenommene Polaritaet des Messdatums falsch,
c     ggf. Korrektur
            if (lpol) then
               imdat  = imdat + dble(idat)*pi
               dat(i) = dcmplx(redat,imdat)

               write(10,'(i4,a19)')
     1              i,' : correct polarity'
               if (.not.lsetup) wdfak(i)=0
            else
               wdfak(i) = 0
            end if

         else if (idat.ne.0.and.isig.ne.0) then

c     Polaritaetswechsel
            vnr(i)    = elec3*10000 + elec4
            imsig     = imsig + dble(isig)*pi
            sigmaa(i) = dcmplx(resig,imsig)
            volt(i)   = cdexp(-sigmaa(i))
            imdat     = imdat + dble(idat)*pi
            dat(i)    = dcmplx(redat,imdat)

            write(10,'(i4,a18)')
     1           i,' : change polarity'

         end if
      end do

c     Ggf. Ausgabe
      do i=1,nanz
         if (wdfak(i).eq.0) write(10,'(i4,a11)') i,' : excluded'
      end do

c     'run.ctr' schliessen
      close(10)

      return
      end

      subroutine bsigm0(kanal,dstart)
      
c     Unterprogramm zum Belegen des Startmodells.

c     Andreas Kemna                                            02-May-1995
c     Letzte Aenderung   15-Jan-2001

c.....................................................................

      USE femmod
      USE datmod
      USE invmod

      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'sigma.fin'
c     diff+<
      INCLUDE 'model.fin'
      INCLUDE 'konv.fin'
c     diff+>

c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Kanalnummer
      integer         * 4     kanal

c     Dateiname
      character       * 80    dstart

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Indexvariable
      integer         * 4     i

c     Hilfsvariablen
      integer         * 4     idat
      real            * 8     dum,dum2

c     Real-, Imaginaerteil
      real            * 8     redat,imdat

c     Pi
      real            * 8     pi

c.....................................................................

      pi = dacos(-1d0)

c     diff-        if (lstart) then
c     diff+<
      if (ldiff) then
         do i=1,elanz
            sigma(i) = cdexp(m0(mnr(i)))
         end do
      else if (lstart) then
c     diff+>
c     'sigma' aus 'dstart' einlesen
         call rsigma(kanal,dstart)
         if (errnr.ne.0) goto 999

      else if (lrho0) then

c     'sigma' gemaess 'bet0', 'pha0' belegen
         do i=1,elanz
            sigma(i) = dcmplx( dcos(pha0/1d3)/bet0 ,
     1           -dsin(pha0/1d3)/bet0 )
         end do

      else if (lbeta) then
         
c     Geometriefaktoren der Messungen bestimmen ("Standard-Geometrie")
         call bkfak()
         if (errnr.ne.0) goto 999

         sigma0 = dcmplx(0d0)
         dum    = 0d0
         
         do i=1,nanz
            
c     Phase lokal korrigieren
c     (entspricht hier "lpol=.true.", aber anders nicht moeglich)
            imdat = dimag(dat(i))
            
            if (imdat.gt.pi/2d0) then
               idat = -1
            else if (imdat.le.-pi/2d0) then
               idat = 1
            else
               idat = 0
            end if

            if (idat.ne.0) imdat=imdat+dble(idat)*pi
            
c     Von "transfer resistance" auf scheinbaren Widerstand umrechnen
            redat = dble(dat(i))-dlog(dabs(kfak(i)))

c     Werte gewichtet mitteln
            dum2   = dsqrt(wmatd(i))*dble(wdfak(i))
            sigma0 = sigma0 + dcmplx(redat,imdat)*dcmplx(dum2)
            dum    = dum + dum2
         end do
         
c     Ggf. Fehlermeldung
         if (dabs(dum).eq.0d0) then
            fetxt = ' '
            errnr = 99
            goto 999
         end if
         
c     'sigma' belegen
         do i=1,elanz
            sigma(i) = cdexp(sigma0/dcmplx(dum))
         end do
         
      else

c     Fehlermeldung
         fetxt = ' '
         errnr = 100
         goto 999
         
      end if
      
c$$$  do i=1,elanz
c$$$  print*,sigma(i),m0(i),cdexp(m0(i))
c$$$  END DO
c     Referenzleitfaehigkeit 'sigma0' bestimmen
      call refsig()

      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 999  return

      end

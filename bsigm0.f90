subroutine bsigm0(kanal,dstart)

!!!$     Unterprogramm zum Belegen des Startmodells.

!!!$     Andreas Kemna                                            02-May-1995
!!!$     Letzte Aenderung   15-Jan-2001

!!!$.....................................................................

  USE femmod
  USE datmod
  USE invmod
  USE sigmamod
  USE modelmod
  USE elemmod
  USE errmod
  USE konvmod

  IMPLICIT none

!!!$     diff+<
!!!$     diff+>

!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Kanalnummer
  INTEGER (KIND = 4)  ::  kanal

!!!$     Dateiname
  CHARACTER (80) ::   dstart

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Indexvariable
  INTEGER (KIND = 4)  ::   i

!!!$     Hilfsvariablen
  INTEGER (KIND = 4)  ::   idat
  REAL (KIND(0D0))    ::   dum,dum2

!!!$     Real-, Imaginaerteil
  REAL (KIND(0D0))    ::   redat,imdat

!!!$     Pi
  REAL (KIND(0D0))    ::   pi

!!!$.....................................................................

  pi = dacos(-1d0)

!!!$     diff-        if (lstart) then
!!!$     diff+<
  if (ldiff) then
     do i=1,elanz
        sigma(i) = cdexp(m0(mnr(i)))
     end do
  else if (lstart) then
!!!$     diff+>
!!!$     'sigma' aus 'dstart' einlesen
     call rsigma(kanal,dstart)
     if (errnr.ne.0) goto 999

  else if (lrho0) then

!!!$     'sigma' gemaess 'bet0', 'pha0' belegen
     do i=1,elanz
        sigma(i) = dcmplx( dcos(pha0/1d3)/bet0 , -dsin(pha0/1d3)/bet0 )
     end do

  else if (lbeta) then

!!!$     Geometriefaktoren der Messungen bestimmen ("Standard-Geometrie")
     call bkfak()
     if (errnr.ne.0) goto 999

     sigma0 = dcmplx(0d0)
     dum    = 0d0

     do i=1,nanz

!!!$     Phase lokal korrigieren
!!!$     (entspricht hier "lpol=.true.", aber anders nicht moeglich)
        imdat = dimag(dat(i))

        if (imdat.gt.pi/2d0) then
           idat = -1
        else if (imdat.le.-pi/2d0) then
           idat = 1
        else
           idat = 0
        end if

        if (idat.ne.0) imdat=imdat+dble(idat)*pi

!!!$     Von "transfer resistance" auf scheinbaren Widerstand umrechnen
        redat = dble(dat(i))-dlog(dabs(kfak(i)))

!!!$     Werte gewichtet mitteln
        dum2   = dsqrt(wmatd(i))*dble(wdfak(i))
        sigma0 = sigma0 + dcmplx(redat,imdat)*dcmplx(dum2)
        dum    = dum + dum2
     end do

!!!$     Ggf. Fehlermeldung
     if (dabs(dum).eq.0d0) then
        fetxt = ' '
        errnr = 99
        goto 999
     end if

!!!$     'sigma' belegen
     do i=1,elanz
        sigma(i) = cdexp(sigma0/dcmplx(dum))
     end do

  else

!!!$     Fehlermeldung
     fetxt = ' '
     errnr = 100
     goto 999

  end if

!!!$  do i=1,elanz
!!!$  print*,sigma(i),m0(i),cdexp(m0(i))
!!!$  END DO
!!!$     Referenzleitfaehigkeit 'sigma0' bestimmen
  call refsig()

  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

999 return

end subroutine bsigm0

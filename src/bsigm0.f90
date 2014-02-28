subroutine bsigm0(kanal,dstart)

!!!$     Unterprogramm zum Belegen des Startmodells.

!!!$     Andreas Kemna                                            02-May-1995
!!!$     Letzte Aenderung   15-Jan-2001

!!!$.....................................................................
use alloci, only:prec
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
  REAL (prec)    ::   dum,dum2

!!!$     Real-, aimaginaerteil
  REAL (prec)    ::   redat,imdat

!!!$     Pi
  REAL (prec)    ::   pi

!!!$.....................................................................

  pi = dacos(-1d0)

!!!$     diff-        if (lstart) then
!!!$     diff+<
  if (ldiff) then ! m0 was set within rall for this case 
!!!!$ (absolute difference inversion)
     do i=1,elanz
        sigma(i) = EXP(m0(mnr(i)))
     end do
  else if (lstart) then
!!!$     diff+>
!!!$     'sigma' aus 'dstart' einlesen
     call rsigma(kanal,dstart)
     if (errnr.ne.0) goto 999

  else if (lrho0) then

!!!$     'sigma' gemaess 'bet0', 'pha0' belegen
     do i=1,elanz
        sigma(i) = CMPLX( COS(pha0/1d3)/bet0 , -SIN(pha0/1d3)/bet0 )
     end do

  else if (lbeta) then

!!!$     Geometriefaktoren der Messungen bestimmen ("Standard-Geometrie")
     call bkfak()
     if (errnr.ne.0) goto 999

     sigma0 = CMPLX(0d0)
     dum    = 0d0

     do i=1,nanz
!!!$     Phase lokal korrigieren
!!!$     (entspricht hier "lpol=.true.", aber anders nicht moeglich)
        imdat = aimag(dat(i))

        if (imdat.gt.pi/2d0) then
           idat = -1
        else if (imdat.le.-pi/2d0) then
           idat = 1
        else
           idat = 0
        end if

        if (idat.ne.0) THEN
           imdat=imdat+REAL(idat)*pi
           PRINT*,'swapping line',idat,i
        END if
        
!!!$     Von "transfer resistance" auf scheinbaren Widerstand umrechnen
        redat = REAL(dat(i))-LOG(ABS(kfak(i)))

!!!$     Werte gewichtet mitteln
        dum2   = SQRT(wmatd(i))*REAL(wdfak(i))
        sigma0 = sigma0 + CMPLX(redat,imdat)*CMPLX(dum2)
        dum    = dum + dum2
!        print*,'Write:',REAL(REAL(sigma0)),REAL(REAL(dat(i))),REAL(redat),REAL(dum2)


     end do

!!!$     Ggf. Fehlermeldung
     if (ABS(dum).eq.0d0) then
        fetxt = 'unable to find starting value sigma0'
        errnr = 99
        goto 999
     end if

!!!$     'sigma' belegen
!     print*,'Write:',sigma0
     do i=1,elanz
        sigma(i) = EXP(sigma0/CMPLX(dum))
     end do

  else

!!!$     Fehlermeldung
     fetxt = ' '
     errnr = 100
     goto 999

  end if


!!$     Referenzleitfaehigkeit 'sigma0' bestimmen

!!!$ >> RM This is now a must have for mixed boundaries
  IF (lbeta) call refsig()
!!!!$ << RM because the sigma0 is needed
  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

999 return

end subroutine bsigm0

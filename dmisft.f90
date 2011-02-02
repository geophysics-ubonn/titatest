subroutine dmisft(lsetup)

!!!$     Unterprogramm zum Bestimmen des Misfits der Daten.

!!!$     Andreas Kemna                                            01-Mar-1995
!!!$     Letzte Aenderung   15-Jan-2001

!!!$.....................................................................

  USE invmod
  USE datmod
  USE errmod
  USE konvmod
  USE pathmod

  IMPLICIT none


!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Hilfsschalter
  LOGICAL ::     lsetup

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsfelder
  REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: psi
  REAL,DIMENSION(:),ALLOCATABLE   :: eps2
  INTEGER(KIND = 4),DIMENSION(:),ALLOCATABLE :: wdlok

!!!$     Hilfsvariablen
  INTEGER (KIND = 4)  ::     i,idum
  COMPLEX (KIND(0D0)) ::    cdum,cdat,csig
  REAL    ::     dum,norm,norm2

!!!$.....................................................................

!!!$     einfach mal oeffnen falls Ausgabe
  errnr = 1
  fetxt = ramd(1:lnramd)//slash(1:1)//'eps.ctr'
  OPEN(fpeps,file=TRIM(fetxt),STATUS='old',POSITION='append',ERR=1000)
  fetxt = ramd(1:lnramd)//slash(1:1)//'run.ctr'
  OPEN(fprun,file=TRIM(fetxt),STATUS='old',POSITION='append',ERR=1000)
  errnr = 4
11 FORMAT (G12.3,2x,F10.6,2x,I3,2X,2G15.7)
  if ((llam.and..not.lstep).or.lsetup) then
     IF (lip) THEN
        write(fpeps,'(/a,t7,I4)',err=1000)'PIT#',it
     ELSE
        write(fpeps,'(/a,t7,I4)',err=1000)'IT#',it
     END IF
     write(fpeps,'(t8,a,t18,a,t27,a,t35,a,t50,a)',err=1000)&
          'eps','psi','pol','Re(f(m))','Im(f(m))'
  end if

!!!$     RMS-WERTE BERECHNEN
  nrmsd  = 0d0
  betrms = 0d0
  pharms = 0d0
  idum   = 0
  !     get memory for wdlok and psi
  ALLOCATE (wdlok(nanz),psi(nanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation psi'
     errnr = 94
     RETURN
  END IF

  do i=1,nanz
     wdlok(i) = 1

!!!$     Phasen lokal korrigieren
     call chkpo2(dat(i),sigmaa(i),cdat,csig,wdlok(i),lpol)

!!!$     ak Ggf. Daten mit Phase betraglich groesser 200 mrad nicht beruecksichtigen
!!!$     ak (Standardabweichung auf 1d4 hochsetzen)
!!!$     ak           if (dabs(1d3*dimag(csig)).gt.200d0) wmatd(i)=1d-8

!!!$     diff-            cdum = cdat - csig
!!!$     diff+<
     if (.not.ldiff) then
        cdum = cdat - csig
     else
        cdum = cdat - csig - (d0(i) - fm0(i))
     end if
!!!$     diff+>

     if (lip) then
        psi(i) = dsqrt(wmatd(i))*dabs(dimag(cdum))
     else
        psi(i) = dsqrt(wmatd(i))*cdabs(cdum)
     end if

!!!$     Ggf. 'eps_i', 'psi_i' und Hilfsfeld ausgeben
     if ((llam.and..not.lstep).or.lsetup) &
          write(fpeps,11,err=1000) real(1d0/dsqrt(wmatd(i))),&
          real(psi(i)),wdlok(i),REAL(csig),AIMAG(csig)

     idum   = idum   + wdlok(i)
     nrmsd  = nrmsd  + psi(i)*psi(i)*dble(wdlok(i))

     betrms = betrms + wmatdr(i)*dble(wdlok(i)) * &
          dble(cdum)*dble(cdum)

     pharms = pharms + wmatdp(i)*dble(wdlok(i)) * &
          dimag(cdum)*dimag(cdum)

  end do

!!!$     Ggf. Fehlermeldung
  if (idum.eq.0) then
     fetxt = ' '
     errnr = 99
     goto 1000
  end if

  npol   = nanz-idum
  rmssum = nrmsd
  nrmsd  = dsqrt(nrmsd /dble(idum))
  betrms = dsqrt(betrms/dble(idum))
  pharms = dsqrt(pharms/dble(idum))

!!!$     Ggf. ROBUST INVERSION (nach Doug' LaBrecque)
  if (lrobust) then

     !     get memory for wdlok
     ALLOCATE (eps2(nanz),stat=errnr)
     IF (errnr /= 0) THEN
        fetxt = 'Error memory allocation eps2'
        errnr = 94
        RETURN
     END IF
!!!$     'estimated weights' und 1-Normen berechnen
     norm  = 0.
     norm2 = 0.

     do i=1,nanz
        dum     = 1./sqrt(wmatd(i))
        eps2(i) = dum*sqrt(psi(i))
        norm    = norm  + REAL(psi(i))*REAL(wdlok(i))
        norm2   = norm2 + REAL(psi(i))*REAL(wdlok(i))*dum/eps2(i)
     end do

!!!$     'estimated weights' normieren
     do i=1,nanz
        eps2(i) = eps2(i) * norm2/norm
     end do

!!!$     Kleinere Standardabweichung ausschliessen und neue 1-Norm berechnen
     norm2 = 0d0

     do i=1,nanz
        dum     = 1./sqrt(wmatd(i))
        eps2(i) = max(dum,eps2(i))
        norm2   = norm2 + REAL(psi(i))*REAL(wdlok(i))*dum/eps2(i)
     end do

     l1rat = norm/norm2

!!!$     Ggf. neue Wichtungsfaktoren belegen
     if (l1rat.gt.l1min) then
        do i=1,nanz
           dum = 1./(eps2(i) * eps2(i))

!!!$     Ausgabe, falls 'eps_neu' > 1.1 * 'eps_alt'
           if (dum.lt.0.83d0*wmatd(i).and. &
                ((llam.and..not.lstep).or.lsetup)) then

              write(fprun,*) &
                   i,' : increase standard deviation *', &
                   eps2(i),dum,psi(i)
           end if

           wmatd(i) = dum
        end do
     end if
     IF (ALLOCATED(eps2)) DEALLOCATE (eps2)
  end if

  IF (ALLOCATED(wdlok)) DEALLOCATE (wdlok,psi)

  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen
1000 close(fpeps)
  close(fprun)
  return

end subroutine dmisft


!!!$*********************************************************************

subroutine chkpo2(dati,sigi,cdat,csig,wdlok,ldum)

!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Eingabe
  COMPLEX (KIND(0D0)) ::    dati,sigi

!!!$     Ausgabe
  COMPLEX (KIND(0D0)) ::    cdat,csig

!!!$     Schalter
  INTEGER (KIND = 4)  ::     wdlok
  LOGICAL ::     ldum

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsvariablen
  INTEGER (KIND = 4)  ::     idat,isig

!!!$     Real-, Imaginaerteile
  REAL (KIND(0D0))    ::     redat,imdat,resig,imsig

!!!$     Pi
  REAL (KIND(0D0))    ::     pi

!!!$.....................................................................

  pi = dacos(-1d0)

!!!$     Logarithmierte Betraege in den Realteilen,
!!!$     Phasen (in rad) in den Imaginaerteilen
  redat = dble(dati)
  imdat = dimag(dati)
  resig = dble(sigi)
  imsig = dimag(sigi)

!!!$     Phasenbereich checken
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

!!!$     Falls ldum=.true., angenommene Polaritaet des Messdatums falsch,
!!!$     ggf. Korrektur; auf jeden Fall Polaritaetswechsel
     imsig = imsig + dble(isig)*pi
     if (.not.ldum) imdat=imdat-dsign(pi,imdat)

     wdlok = 0

  else if (idat.ne.0.and.isig.eq.0) then

!!!$     Falls ldum=.true., angenommene Polaritaet des Messdatums falsch,
!!!$     ggf. Korrektur
     if (ldum) imdat=imdat+dble(idat)*pi

     wdlok = 0

  else if (idat.ne.0.and.isig.ne.0) then

!!!$     Polaritaetswechsel
     imsig = imsig + dble(isig)*pi
     imdat = imdat + dble(idat)*pi

  end if

!!!$     'cdat' und 'csig' speichern
  cdat = dcmplx(redat,imdat)
  csig = dcmplx(resig,imsig)

  return
end subroutine chkpo2

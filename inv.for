      program inv

c     Hauptprogramm zur Complex-Resistivity-2.5D-Inversion.

c     Belegte Kanaele:  9 - error.dat
c     10 - run.ctr
c     11 - in-/output
c     12 - crtomo.cfg
c     13 - inv.ctr
c     14 - cjg.ctr
c     15 - eps.ctr

c     Andreas Kemna                                            02-May-1995
c     Letzte Aenderung   22-Jul-2007

c.....................................................................

      USE alloci
c     USE portlib

      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'
      INCLUDE 'invhp.fin'
      INCLUDE 'path.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'electr.fin'
      INCLUDE 'waven.fin'
      INCLUDE 'sigma.fin'
      INCLUDE 'dat.fin'
      INCLUDE 'model.fin'
      INCLUDE 'fem.fin'
      INCLUDE 'inv.fin'
      INCLUDE 'konv.fin'

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     SETUP UND INPUT
c     'crtomo.cfg' oeffnen
      fetxt = 'crtomo.cfg'
      errnr = 1
      open(12,file=fetxt,status='old',err=999)

c     Allgemeine Parameter setzen
      kanal  = 11
 5    errnr2 = 0

c     Benoetigte Variablen einlesen
      call rall(kanal,delem,delectr,dstrom,drandb,
c     diff-     1            dsigma,dvolt,dsens,dstart,lsens,lagain)
c     diff+<
     1     dsigma,dvolt,dsens,dstart,dd0,dm0,dfm0,
     1     lsens,lagain)
c     diff+>
      if (errnr.ne.0) goto 999

      PRINT*,'ltri::',ltri

c     Element- und Randelementbeitraege sowie ggf. Konfigurationsfaktoren
c     zur Berechnung der gemischten Randbedingung bestimmen
      call precal()
      if (errnr.ne.0) goto 999

      if (.not.lbeta) then
         lsr = .false.

c     Ggf. Fehlermeldungen
         if (.not.lsink) then
            fetxt = ' '
            errnr = 102
         end if

         if (swrtr.eq.0.and..not.lrandb2) then
            fetxt = ' '
            errnr = 103
         end if

         if (errnr.ne.0) goto 999
      else

c     Ggf. Fehlermeldungen
         if (swrtr.eq.0) then
            fetxt = ' '
            errnr = 106
         end if
      end if

c     Startmodell belegen
      call bsigm0(kanal,dstart)
      if (errnr.ne.0) goto 999        

c     Startparameter setzen
      it     = 0
      itr    = 0
      rmsalt = 0d0
      lsetup = .true.
      lsetip = .false.
      lip    = .false.
      llam   = .false.
      ldlami = .true.
      lstep  = .false.
      lfstep = .false.
      step   = 1d0
      stpalt = 1d0
      alam   = 0d0

c     Kontrolldateien initialisieren
c     diff-        call kont1(delem,delectr,dstrom,drandb)
c     diff+<
      call kont1(delem,delectr,dstrom,drandb,dd0,dm0,dfm0)
c     diff+>
      if (errnr.ne.0) goto 999

c     'sens' zuweisen
      if (ldc) then
         ALLOCATE(sensdc(nanz,manz),stat=errnr)
      else
         ALLOCATE(sens(nanz,manz),stat=errnr)
      end if
      if (errnr.ne.0) then
         errnr = 97 
         goto 999
      end if

c.................................................

c     MODELLING
c     'a', 'hpot' und 'kpot' zuweisen
 10   if (ldc) then
         ALLOCATE(adc((mb+1)*sanz),hpotdc(sanz,eanz),
     1        kpotdc(sanz,eanz,kwnanz),stat=errnr)
      else
         ALLOCATE(a((mb+1)*sanz),hpot(sanz,eanz),
     1        kpot(sanz,eanz,kwnanz),stat=errnr)
      end if
      if (errnr.ne.0) then
         errnr = 97 
         goto 999
      end if

c     Kontrollausgaben
      write(*,'(a11,i3,a2,i3,a25)')
     1     ' Iteration ',it,', ',itr,
     1     ' : Calculating Potentials'

      errnr = 1
      fetxt = ramd(1:lnramd)//slash(1:1)//'run.ctr'
      open(10,file=fetxt,status='unknown',err=999)
 11   read(10,*,end=12)
      goto 11
 12   backspace(10)
      errnr = 4
      write(10,'(a11,i3,a2,i3,a25)',err=999)
     1     ' Iteration ',it,', ',itr,
     1     ' : Calculating Potentials'
      close(10)

      if (ldc) then

c     DC CASE
         do k=1,kwnanz
            do l=1,eanz
               if (lsr.or.lbeta.or.l.eq.1) then

c     Ggf. Potentialwerte fuer homogenen Fall analytisch berechnen
                  if (lsr) call potana(l,k)

c     Kompilation des Gleichungssystems (fuer Einheitsstrom !)
                  call kompadc(l,k)
                  if (errnr.ne.0) goto 999

c     Ggf. Randbedingung beruecksichtigen
                  if (lrandb) call randdc()
                  if (lrandb2) call randb2()

c     Gleichungssystem skalieren
                  call scaldc()
                  if (errnr.ne.0) goto 999

c     Cholesky-Zerlegung der Matrix
                  call choldc()
                  if (errnr.ne.0) goto 999
               else

c     Stromvektor modifizieren
                  call kompbdc(l)
               end if

c     Gleichungssystem loesen
               call vredc()

c     Potentialwerte zurueckskalieren und umspeichern sowie ggf.
c     analytische Loesung addieren
               do j=1,sanz
                  kpotdc(j,l,k) = dble(pot(j)) * fak(j)
                  if (lsr) kpotdc(j,l,k) = kpotdc(j,l,k) +
     1                 dble(pota(j))
                  if (swrtr.eq.0) hpotdc(j,l) = kpotdc(j,l,k)
               end do
            end do
         end do

      else

c     COMPLEX CASE
         do k=1,kwnanz
            do l=1,eanz
               if (lsr.or.lbeta.or.l.eq.1) then

c     Ggf. Potentialwerte fuer homogenen Fall analytisch berechnen
                  if (lsr) call potana(l,k)

c     Kompilation des Gleichungssystems (fuer Einheitsstrom !)
                  call kompab(l,k)
                  if (errnr.ne.0) goto 999

c     Ggf. Randbedingung beruecksichtigen
                  if (lrandb) call randb()
                  if (lrandb2) call randb2()

c     Gleichungssystem skalieren
                  call scalab()
                  if (errnr.ne.0) goto 999

c     Cholesky-Zerlegung der Matrix
                  call chol()
                  if (errnr.ne.0) goto 999
               else

c     Stromvektor modifizieren
                  call kompb(l)
               end if

c     Gleichungssystem loesen
               call vre()

c     Potentialwerte zurueckskalieren und umspeichern sowie ggf.
c     analytische Loesung addieren
               do j=1,sanz
                  kpot(j,l,k) = pot(j) * dcmplx(fak(j))
                  if (lsr) kpot(j,l,k) = kpot(j,l,k) + pota(j)
                  if (swrtr.eq.0) hpot(j,l) = kpot(j,l,k)
               end do
            end do
         end do

      end if

c     Ggf. Ruecktransformation der Potentialwerte
      if (swrtr.eq.1) call rtrafo()

c     Spannungswerte berechnen
      call bvolti()
      if (errnr.ne.0) goto 999

c     'a' und 'hpot' freigeben
      if (ldc) then
         DEALLOCATE(adc,hpotdc)
      else
         DEALLOCATE(a,hpot)
      end if

      if (lsetup) then

c     Ggf. background auf ratio-Daten "multiplizieren"
         if (lratio) then
            do j=1,nanz
               dat(j) = dat(j) + sigmaa(j)
            end do
         end if

c     Polaritaeten checken
         call chkpol(lsetup)
      end if

c     Daten-RMS berechnen
      call dmisft(lsetup)
      if (errnr.ne.0) goto 999

c     'nrmsd=0' ausschliessen
      if (nrmsd.lt.1d-12) nrmsd=nrmsdm*(1d0-mqrms)

c     tst
c     tst        if (lfphai) then
c     tst            llam = .true.
c     tst            if (.not.lip) nrmsd = 1d0
c     tst        end if

c.............................

c     Kontrollvariablen ausgeben
      call kont2(lsetup)
      if (errnr.ne.0) goto 999

c     ABBRUCHBEDINGUNGEN
      if (llam.and..not.lstep) then

c     Polaritaeten checken
         call chkpol(lsetup)

c     Wiederholt minimale step-length ?
         if (stpalt.eq.0d0) errnr2=92

c     Keine Verbesserung des Daten-RMS ?
         if (dabs(1d0-rmsalt/nrmsd).le.mqrms) errnr2=81

c     Minimaler Daten-RMS erreicht ?
c     tst            if (dabs(1d0-nrmsd/nrmsdm).le.mqrms) errnr2=80
         if (dabs(1d0-nrmsd/nrmsdm).le.mqrms.and.ldlamf) errnr2=80

c     Maximale Anzahl an Iterationen ?
         if (it.ge.itmax) errnr2=79

c     Ggf. abbrechen oder "final phase improvement"
         if (errnr2.ne.0) then
            if (lfphai.and.errnr2.ne.79) then
               errnr2 = 0
c     ak
c     Widerstandsverteilung und modellierte Daten ausgeben
               call wout(kanal,dsigma,dvolt)
               if (errnr.ne.0) goto 999

c     Kontrollausgaben
               write(*,'(a24)') ' Final phase improvement'

               errnr = 1
               fetxt = ramd(1:lnramd)//slash(1:1)//'run.ctr'
               open(10,file=fetxt,status='old',err=999)
 13            read(10,*,end=14)
               goto 13
 14            backspace(10)
               errnr = 4
               write(10,'(a24)',err=999)
     1              ' Final phase improvement'
               close(10)

               errnr = 1
               fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
               open(13,file=fetxt,status='old',err=999)
 15            read(13,*,end=16)
               goto 15
 16            backspace(13)
               errnr = 4
               write(13,'(a48,a48,a12)',err=999)
     1              '------------------------------------------------',
     1              '------------------------------------------------',
     1              '------------'
               write(13,*,err=999)
               close(13)

c     Wichtungsfeld umspeichern
               do j=1,nanz
                  wmatd(j) = wmatdp(j)
               end do

               lip    = .true.
               lsetip = .true.
               lfphai = .false.
               llam   = .false.
               ldlami = .true.
               lfstep = .true.
               step   = 1d0

c     ak
               do j=1,elanz
                  sigma(j) = dcmplx(
     1                 dcos(pha0/1d3)*cdabs(sigma(j)) ,
     1                 -dsin(pha0/1d3)*cdabs(sigma(j)) )
               end do
c     ak

c     Daten-RMS berechnen
               call dmisft(lsetip)
               if (errnr.ne.0) goto 999

c     Kontrollvariablen ausgeben
               call kont2(lsetip)
               if (errnr.ne.0) goto 999
            else
               goto 30
            end if
         else
c     ak
c     Widerstandsverteilung und modellierte Daten ausgeben
            call wout(kanal,dsigma,dvolt)
            if (errnr.ne.0) goto 999
         end if
      end if

      if ((llam.and..not.lstep).or.lsetup.or.lsetip) then

c     Iterationsindex hochzaehlen
         it = it+1

c     Parameter zuruecksetzen
         itr    = 0
         rmsreg = 0d0
         dlam   = 1d0
         dlalt  = 1d0
         ldlamf = .false.

c     Daten-RMS speichern
         rmsalt = nrmsd

c     Felder speichern
         do j=1,elanz
            sigma2(j) = sigma(j)
         end do
         do j=1,nanz
            if (lrobust) wmatd2(j)=wmatd(j)
            sgmaa2(j) = sigmaa(j)
         end do

c     Kontrollausgaben
         write(*,'(a11,i3,a2,i3,a28)')
     1        ' Iteration ',it,', ',itr,
     1        ' : Calculating Sensitivities'

         errnr = 1
         fetxt = ramd(1:lnramd)//slash(1:1)//'run.ctr'
         open(10,file=fetxt,status='old',err=999)
 17      read(10,*,end=18)
         goto 17
 18      backspace(10)
         errnr = 4
         write(10,'(a11,i3,a2,i3,a28)')
     1        ' Iteration ',it,', ',itr,
     1        ' : Calculating Sensitivities'
         close(10)

c     SENSITIVITAETEN berechnen
         if (ldc) then
            call bsendc()
         else
            call bsensi()
         end if

         if (lsetup) then
c     akc Ggf. Summe der Sensitivitaeten aller Messungen ausgeben
c     if (lsens) then
c     if (ldc) then
c     call bbsedc(kanal,dsens)
c     else
c     call bbsens(kanal,dsens)
c     end if
c     if (errnr.ne.0) goto 999
c     end if
c     ak
c     Rauhigkeitsmatrix belegen
            PRINT*,'first iter..',ltri
            IF (ltri==0) THEN
               call bsmatm()
            ELSE IF (ltri==1) THEN
               PRINT*,'Here:: Triangulation regularization'
               CALL bsmatmtri
            END IF 
c     ak



c     open(21,file='smatm.dat',
c     1                  access='sequential',form='unformatted')
c     do j=1,mmax
c     c                    write(21) (smatm(j,k),k=1,3)
c     read(21) (smatm(j,k),k=1,3)
c     end do
c     close(21)
c     ak
         end if
      else
c     Felder zuruecksetzen
         do j=1,elanz
            sigma(j)  = sigma2(j)
         end do
         do j=1,nanz
            if (lrobust) wmatd(j)=wmatd2(j)
            sigmaa(j) = sgmaa2(j)
         end do
      end if

c     'kpot' freigeben
      if (ldc) then
         DEALLOCATE(kpotdc)
      else
         DEALLOCATE(kpot)
      end if

c     REGULARISIERUNG / STEP-LENGTH einstellen
      if (.not.lstep) then
         if (llam) then

c     "Regularisierungsschleife" initialisieren und step-length zuruecksetzen
            llam = .false.
            step = 1d0
         else

c     Regularisierungsindex hochzaehlen
            itr = itr+1

            if ((((nrmsd.lt.rmsreg.and.itr.le.nlam).or.
     1           (dlam.gt.1d0.and.itr.le.nlam)).and.
     1           (.not.ldlamf.or.dlalt.le.1d0).and.
     1           dabs(1d0-rmsreg/nrmsdm).gt.mqrms).or.
     1           rmsreg.eq.0d0) then

c     Regularisierungsparameter bestimmen
               if (lsetup.or.lsetip) then

c     Kontrollausgabe
                  write(*,'(a11,i3,a2,i3,a43)')
     1                 ' Iteration ',it,', ',itr,
     1                 ' : Calculating 1st regularization parameter'
                  call blam0()
                  lam = lammax
c     ak Model EGS2003, ERT2003                        call blam0()
c     ak Model EGS2003, ERT2003                        lam = lammax
c     ak                        lam = 1d4
               else
                  dlalt = dlam
                  if (ldlami) then
                     ldlami = .false.
                     alam   = dmax1(dabs(dlog(nrmsd/nrmsdm)),
     1                    dlog(1d0+mqrms))
                     dlam   = fstart
                  else
                     alam = dmax1(alam,dabs(dlog(nrmsd/nrmsdm)))
                     dlam = dlog(fstop)*
     1                    sign(1d0,dlog(nrmsd/nrmsdm))+
     1                    dlog(fstart/fstop)*
     1                    dlog(nrmsd/nrmsdm)/alam
                     dlam = dexp(dlam)
                  end if
                  lam = lam*dlam
                  if (dlalt.gt.1d0.and.dlam.lt.1d0) ldlamf=.true.
c     tst                        if (dlam.gt.1d0) lfstep=.true.
c     ak Model EGS2003
                  if (dlam.gt.1d0) lrobust=.false.
               end if
            else

c     Regularisierungsparameter zuruecksetzen und step-length verkleinern
               llam = .true.
               lam  = lam/dlam

               if (lfstep) then
                  lfstep = .false.
               else
                  lstep = .true.
                  step  = 5d-1
               end if
            end if

c     Ggf. Daten-RMS speichern
            if (lsetup.or.lsetip) then
               lsetup = .false.
               lsetip = .false.
            else
               if (.not.lstep) rmsreg=nrmsd
            end if
         end if
      else
         lstep = .false.

c     Parabolische Interpolation zur Bestimmung der optimalen step-length
         call parfit(rmsalt,nrmsd,rmsreg,nrmsdm,stpmin)

         if (step.eq.stpmin.and.stpalt.eq.stpmin) then

c     Nach naechstem Modelling abbrechen
            stpalt = 0d0
         else

c     Step-length speichern
            stpalt = step
         end if
      end if

c     Kontrollausgaben
      write(*,'(a11,i3,a2,i3,a11)')
     1     ' Iteration ',it,', ',itr,
     1     ' : Updating'

      errnr = 1
      fetxt = ramd(1:lnramd)//slash(1:1)//'run.ctr'
      open(10,file=fetxt,status='old',err=999)
 19   read(10,*,end=20)
      goto 19
 20   backspace(10)
      errnr = 4
      write(10,'(a11,i3,a2,i3,a11)')
     1     ' Iteration ',it,', ',itr,
     1     ' : Updating'
      close(10)

c     UPDATE anbringen
      call update(dpar2,cgres2)

c     Leitfaehigkeiten belegen und Roughness bestimmen
      IF (ltri==0) THEN
         PRINT*,'Old roughness'
         call brough()
      ELSE IF (ltri==1) THEN
         PRINT*,'Triangulation roughness'         
         CALL broughtri
      END IF

c     Ggf. Referenzleitfaehigkeit bestimmen
      if (lsr) call refsig()

c     Neues Modelling
      goto 10

c.................................................

c     OUTPUT
 30   call wout(kanal,dsigma,dvolt)
      if (errnr.ne.0) goto 999

c     Ggf. Summe der Sensitivitaeten aller Messungen ausgeben
      if (lsens) then
         if (ldc) then
            call bbsedc(kanal,dsens)
         else
            call bbsens(kanal,dsens)
         end if
         if (errnr.ne.0) goto 999
      end if

c     Kontrollausgaben
      errnr = 1
      fetxt = ramd(1:lnramd)//slash(1:1)//'run.ctr'
      open(10,file=fetxt,status='old',err=999)
 31   read(10,*,end=32)
      goto 31
 32   backspace(10)
      errnr = 4

      if (errnr2.eq.92) then
         write(*,'(a22,a31)') ' Iteration terminated:',
     1        ' Min. step-length for 2nd time.'

         write(10,'(a22,a31)',err=999) ' Iteration terminated:',
     1        ' Min. step-length for 2nd time.'
      else if (errnr2.eq.80) then
         write(*,'(a22,a10)') ' Iteration terminated:',
     1        ' Min. RMS.'

         write(10,'(a22,a10)',err=999) ' Iteration terminated:',
     1        ' Min. RMS.'
      else if (errnr2.eq.81) then
         write(*,'(a22,a24)') ' Iteration terminated:',
     1        ' Min. rel. RMS decrease.'

         write(10,'(a22,a24)',err=999) ' Iteration terminated:',
     1        ' Min. rel. RMS decrease.'
      else if (errnr2.eq.79) then
         write(*,'(a22,a19)') ' Iteration terminated:',
     1        ' Max. # iterations.'

         write(10,'(a22,a19)',err=999) ' Iteration terminated:',
     1        ' Max. # iterations.'
      end if

c     Run-time abfragen und ausgeben
      izeit     = etime(tazeit)
      izeit     = izeit/60.
      tazeit(2) = amod(izeit,60.)
      tazeit(1) = (izeit-tazeit(2))/60.

      write(*,'(a10,i3,a3,i3,a3)')
     1     ' CPU time:',int(tazeit(1)),'hrs',int(tazeit(2)),'min'
      write(10,'(a10,i3,a3,i3,a3)',err=999)
     1     ' CPU time:',int(tazeit(1)),'hrs',int(tazeit(2)),'min'
      close(10)

c     Kontrolldateien schliessen
      close(14)
      close(15)

c     'sens' und 'kpot' freigeben
      if (ldc) then
         DEALLOCATE(sensdc,kpotdc)
      else
         DEALLOCATE(sens,kpot)
      end if

c     Ggf. weiteren Datensatz invertieren
      if (lagain) goto 5

c     'crtomo.cfg' schliessen
      close(12)

      stop ' '

c.....................................................................

c     Fehlermeldung
 999  open(9,file='error.dat',status='replace')
      errflag = 2
      write(9,'(a80,i3,i1)') fetxt,errnr,errflag
      close(9)
      stop ' '

      end

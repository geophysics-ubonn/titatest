      program inv

c     Hauptprogramm zur Complex-Resistivity-2.5D-Inversion.

c     Belegte Kanaele:  
c     9 - error.dat   -> fperr
c     10 - run.ctr    -> fprun
c     11 - in-/output -> kanal
c     12 - crtomo.cfg -> fpcfg
c     13 - inv.ctr    -> fpinv
c     14 - cjg.ctr    -> fpcjg
c     15 - eps.ctr    -> fpeps
c
c     Andreas Kemna                                        02-May-1995
c     Letzte Aenderung                                     Jul-2010

c.....................................................................

      USE alloci
      USE tic_toc
      USE femmod
      USE datmod
      USE invmod
      USE cjgmod
      USE sigmamod
      USE electrmod
      USE modelmod
      USE elemmod
      USE wavenmod
      USE randbmod
      USE konvmod
      USE errmod
      USE pathmod
      USE bsmatm_mod
      USE bmcm_mod
      USE brough_mod

c     USE portlib

      IMPLICIT none

      INCLUDE 'invhp.fin'

      CHARACTER(256)         :: ftext
      INTEGER                :: c1,i
      REAL(KIND(0D0))        :: lamalt
      LOGICAL                :: ols
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     SETUP UND INPUT
c     'crtomo.cfg' oeffnen
      fetxt = 'crtomo.cfg'
      errnr = 1
c     Kanal nummern belegen.. die variablen sind global!
      fperr = 9 
      fprun = 10
      kanal = 11 
      fpcfg = 12 
      fpinv = 13 
      fpcjg = 14 
      fpeps = 15 

      open(fpcfg,file=fetxt,status='old',err=999)
      
c     Allgemeine Parameter setzen
 5    errnr2 = 0
c     Benoetigte Variablen einlesen
      call rall(kanal,delem,delectr,dstrom,drandb,
c     diff-     1            dsigma,dvolt,dsens,dstart,lsens,lagain)
c     diff+<
     1     dsigma,dvolt,dsens,dstart,dd0,dm0,dfm0,lagain)
c     diff+>
      if (errnr.ne.0) goto 999


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
      ALLOCATE (sigma(elanz),sigma2(elanz),stat=errnr)
      IF (errnr /= 0) THEN
         fetxt = 'Error memory allocation fem sigma'
         errnr = 94
         goto 999
      END IF
      call bsigm0(kanal,dstart)
      if (errnr.ne.0) goto 999

c     get memory for model parameters
      ALLOCATE (par(manz),dpar(manz),dpar2(manz),stat=errnr)
      IF (errnr /= 0) THEN
         fetxt = 'Error memory allocation model/update data'
         errnr = 94
         goto 999
      END IF
c     get memory for CG data storage of residuums
      ALLOCATE (cgres(ncgmax+2),cgres2(ncgmax+2),stat=errnr)
      IF (errnr /= 0) THEN
         fetxt = 'Error memory allocation cgres data'
         errnr = 94
         goto 999
      END IF

c     Startparameter setzen
      it     = 0
      itr    = 0
      rmsalt = 0d0
      lamalt = 1d0
      bdpar = 1d0
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

c     Kontrolldateien oeffnen
      errnr = 1

      fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
      open(fpinv,file=fetxt,status='replace',err=999)
      close(fpinv)
      fetxt = ramd(1:lnramd)//slash(1:1)//'run.ctr'
      open(fprun,file=fetxt,status='replace',err=999)
c$$$  close(fprun) muss geoeffnet bleiben da sie staendig beschrieben wird
      fetxt = ramd(1:lnramd)//slash(1:1)//'cjg.ctr'
      open(fpcjg,file=fetxt,status='replace',err=999)
      close(fpcjg)
      fetxt = ramd(1:lnramd)//slash(1:1)//'eps.ctr'
      open(fpeps,file=fetxt,status='replace',err=999)
      IF (ldc) THEN
         WRITE (fpeps,'(a)')'1/eps_r'
         WRITE (fpeps,'(G10.3)')(sqrt(wmatdr(i)),i=1,nanz)
      ELSE
         WRITE (fpeps,'(3(a,10x))')'1/eps','1/eps_r','1/eps_p'
         WRITE (fpeps,'(3(G12.5,2x))')
     1        (sqrt(wmatd(i)),sqrt(wmatdr(i)),sqrt(wmatdp(i)),i=1,nanz)
      END IF
      close(fpeps)
      errnr = 4

c     Kontrolldateien initialisieren
c     diff-        call kont1(delem,delectr,dstrom,drandb)
c     diff+<
      call kont1(delem,delectr,dstrom,drandb,dd0,dm0,dfm0)
c     diff+>
      if (errnr.ne.0) goto 999
!!$      CALL SYSTEM('sleep 1000')
c     'sens' zuweisen
      if (ldc) then
         ALLOCATE(sensdc(nanz,manz),kpotdc(sanz,eanz,kwnanz),stat=errnr)
      else
         ALLOCATE(sens(nanz,manz),kpot(sanz,eanz,kwnanz),stat=errnr)
      end if
      if (errnr.ne.0) then
         fetxt = 'allocation problem sens and kpot'
         errnr = 97 
         goto 999
      end if

c-------------
c     get current time
      CALL tic(c1)
c.................................................

c     MODELLING
c     'a', 'hpot' und 'kpot' zuweisen
 10   if (ldc) then
         ALLOCATE(adc((mb+1)*sanz),hpotdc(sanz,eanz),bdc(sanz),
     1        stat=errnr)
      else
         ALLOCATE(a((mb+1)*sanz),hpot(sanz,eanz),b(sanz),stat=errnr)
      end if
      if (errnr.ne.0) then
         fetxt = 'allocation problem a and hpot'
         errnr = 97 
         goto 999
      end if
      IF (.NOT.ALLOCATED (pot)) THEN
         ALLOCATE(pot(sanz),pota(sanz),fak(sanz),stat=errnr)
         if (errnr.ne.0) then
            fetxt = 'allocation problem pot to fak'
            errnr = 97 
            goto 999
         end if
      END IF

c     Kontrollausgaben
      WRITE (*,'(a60)',ADVANCE='no')ACHAR(13)//''
      write(*,'(a,i3,a,i3,a)',ADVANCE='no')
     1     ACHAR(13)//' Iteration ',it,', ',itr,
     1     ' : Calculating Potentials'


      write(fprun,'(a,i3,a,i3,a)')
     1     ' Iteration ',it,', ',itr,
     1     ' : Calculating Potentials'

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
         DEALLOCATE(adc,hpotdc,bdc)
      else
         DEALLOCATE(a,hpot,b)
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
         lamalt = lam 
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

c     Minimal stepsize erreicht ?
         IF (bdpar < bdmin) THEN
            errnr2=109
            WRITE (ftext,*)'check stepsize',bdpar,it,itr
         END IF

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

               write(fprun,'(a24)',err=999)
     1              ' Final phase improvement'

               fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
               open(fpinv,file=fetxt,status='old',
     1              POSITION='append',err=999)
               write(fpinv,'(/a/)',err=999)
     1              '------------------------------------------------'//
     1              '------------------------------------------------'//
     1              '-----------------'
               close(fpinv)

c     Wichtungsfeld umspeichern
               wmatd(1:nanz) = wmatdp(1:nanz)

               lip    = .true.
               lsetip = .true. ! 
               lfphai = .false.
               llam   = .false.
               ldlami = .true.
               lfstep = .true.
               step   = 1d0

c     ak
               fetxt = 'cp -f inv.lastmod inv.lastmod_rho'
               CALL SYSTEM (TRIM(fetxt))
               IF (lffhom) THEN
                  write(*,*)
     1                 ' ******* Restarting phase model ********'
                  write(fprun,*)
     1                 ' ******* Restarting phase model ********'
                  do j=1,elanz
                     sigma(j) = dcmplx(
     1                    dcos(pha0/1d3)*cdabs(sigma(j)) ,
     1                    -dsin(pha0/1d3)*cdabs(sigma(j)) )
                  end do
                  
                  GOTO 10       ! neues calc
               END IF
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

c     Lambda speichen
         IF (it>1) lamalt = lam ! mindestens einmal konvergiert

c     Felder speichern
         sigma2 = sigma

         sgmaa2 = sigmaa

         IF (lrobust) wmatd2 = wmatd

c     Kontrollausgaben
         write (*,'(a,i3,a,i3,a)',ADVANCE='no')
     1        ACHAR(13)//' Iteration ',it,', ',itr,
     1        ' : Calculating Sensitivities'
         
         write(fprun,'(a,i3,a,i3,a)')' Iteration ',it,', ',itr,
     1        ' : Calculating Sensitivities'

c     SENSITIVITAETEN berechnen
         if (ldc) then
            call bsendc()
         else
            call bsensi()
         end if
         
c     evtl   Rauhigkeitsmatrix belegen
         IF (lsetup.OR.(ltri>3.AND.ltri<10)) CALL bsmatm(it) 
                        
      else
c     Felder zuruecksetzen
         sigma = sigma2

         sigmaa = sgmaa2

         IF (lrobust) wmatd = wmatd2

      end if
      IF (itmax == 0) THEN ! only precalcs..
         errnr2 = 109
         print*,'Only precalcs'
         GOTO 30
      END IF 
c$$$c     'kpot' freigeben wird nicht mehr freigegeben
c$$$      if (ldc) then
c$$$         DEALLOCATE(kpotdc)
c$$$      else
c$$$         DEALLOCATE(kpot)
c$$$      end if

c     REGULARISIERUNG / STEP-LENGTH einstellen
      if (.not.lstep) then
         if (llam) then

c     "Regularisierungsschleife" initialisieren und step-length zuruecksetzen
            llam = .false.
            step = 1d0
         else

c     Regularisierungsindex hochzaehlen
            itr = itr+1
            if (((((nrmsd.lt.rmsreg.and.itr.le.nlam).or.
     1           (dlam.gt.1d0.and.itr.le.nlam)).and.
     1           (.not.ldlamf.or.dlalt.le.1d0).and.
     1           dabs(1d0-rmsreg/nrmsdm).gt.mqrms).or.
     1           (rmsreg.eq.0d0)).AND.
     1           (bdpar >= bdmin)) then
c     Regularisierungsparameter bestimmen
               if (lsetup.or.lsetip) then

c     Kontrollausgabe
                  IF (llamf) THEN
                     lam = lamfix
                  ELSE
                     write(*,'(a,i3,a,i3,a)',ADVANCE='no')
     1                    ACHAR(13)//' Iteration ',it,', ',itr,
     1                    ' : Calculating 1st regularization parameter'
                     write(fprun,'(a,i3,a,i3,a)',ADVANCE='no')
     1                    ' Iteration ',it,', ',itr,
     1                    ' : Calculating 1st regularization parameter'
                     call blam0()
                     WRITE (*,'(a,G10.2)',ADVANCE='no')'lam_0:: ',lammax
                     WRITE (fprun,'(a,G10.2)')'lam_0 ',lammax
                     lam = lammax
c     ak Model EGS2003, ERT2003                        call blam0()
c     ak Model EGS2003, ERT2003                        lam = lammax
c     ak                        lam = 1d4
                  END IF
               else
                  IF (llamf) THEN
                     lam = lamfix
                  ELSE
                     dlalt = dlam
                     if (ldlami) then
                        ldlami = .false.
                        alam   = dmax1(dabs(dlog(nrmsd/nrmsdm)),
     1                       dlog(1d0+mqrms))
                        dlam   = fstart
                     else
                        alam = dmax1(alam,dabs(dlog(nrmsd/nrmsdm)))
                        dlam = dlog(fstop)*
     1                       sign(1d0,dlog(nrmsd/nrmsdm))+
     1                       dlog(fstart/fstop)*
     1                       dlog(nrmsd/nrmsdm)/alam
                        dlam = dexp(dlam)
                     end if
                     lam = lam*dlam
                     if (dlalt.gt.1d0.and.dlam.lt.1d0) ldlamf=.true.
c     tst                        if (dlam.gt.1d0) lfstep=.true.
c     ak Model EGS2003
                     if (dlam.gt.1d0) lrobust=.false.
                  END IF
               end if
            else

c     Regularisierungsparameter zuruecksetzen und step-length verkleinern
               llam = .true.
               IF (llamf) THEN
                  lam = lamfix
               ELSE
                  lam  = lam/dlam
               END IF
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

         if (step.eq.stpmin.and.stpalt.eq.stpmin)then

c     Nach naechstem Modelling abbrechen
            stpalt = 0d0
         else

c     Step-length speichern
            stpalt = step
         end if
      end if
c     Kontrollausgaben
      write(*,'(a,i3,a,i3,a)',ADVANCE='no')
     1     ACHAR(13)//' Iteration ',it,', ',itr,
     1     ' : Updating'

      write(fprun,*)' Iteration ',it,', ',itr,
     1     ' : Updating'

c Modell parameter mit aktuellen Leitfaehigkeiten belegen

      CALL bpar
      if (errnr.ne.0) goto 999

c     UPDATE anbringen
      call update
      fetxt =''
      if (errnr.ne.0) goto 999

c Leitfaehigkeiten mit verbessertem Modell belegen
      CALL bsigma
      if (errnr.ne.0) goto 999

c     Roughness bestimmen
      CALL brough
      
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
      IF (lvario) CALL bvariogram ! calculate experimental variogram
      
      IF (lcov1) CALL buncert (kanal,lamalt)

c     Kontrollausgaben

      if (errnr2.eq.92) then
         
         write(*,'(a22,a31)') ' Iteration terminated:',
     1        ' Min. step-length for 2nd time.'
         
         write(fprun,'(a22,a31)',err=999) ' Iteration terminated:',
     1        ' Min. step-length for 2nd time.'
      else if (errnr2.eq.80) then
         write(*,'(a22,a10)') ' Iteration terminated:',
     1        ' Min. RMS.'
         
         write(fprun,'(a22,a10)',err=999) ' Iteration terminated:',
     1        ' Min. RMS.'
      else if (errnr2.eq.81) then
         write(*,'(a22,a24)') ' Iteration terminated:',
     1        ' Min. rel. RMS decrease.'
         
         write(fprun,'(a22,a24)',err=999) ' Iteration terminated:',
     1        ' Min. rel. RMS decrease.'
      else if (errnr2.eq.79) then
         write(*,'(a22,a19)') ' Iteration terminated:',
     1        ' Max. # iterations.'
         
         write(fprun,'(a22,a19)',err=999) ' Iteration terminated:',
     1        ' Max. # iterations.'
         
      else if (errnr2.eq.109) then
         write(*,'(a)') ' Iteration terminated:'//
     1        ' Min. model changes reached'
         
         write(fprun,'(a)',err=999) ' Iteration terminated:'//
     1        ' Min. model changes reached'
         
      end if
      
c     Run-time abfragen und ausgeben
      fetxt = ' CPU time: '
      CALL toc(c1,fetxt)
      write(fprun,'(a)',err=999)TRIM(fetxt)
      
c     Kontrolldateien schliessen
      close(fprun)
      close(fpinv)
      close(fpcjg)
      close(fpeps)
      
c     'sens' und 'kpot' freigeben
      if (ldc) then
         DEALLOCATE(sensdc,kpotdc)
      else
         DEALLOCATE(sens,kpot)
      end if
      IF (ALLOCATED (smatm)) DEALLOCATE (smatm)
      IF (ALLOCATED (pot)) DEALLOCATE (pot,pota,fak)

      IF (ALLOCATED (snr)) DEALLOCATE (snr,sx,sy)
      IF (ALLOCATED (typ)) DEALLOCATE (typ,nelanz,selanz)
      IF (ALLOCATED (nrel)) DEALLOCATE (nrel,rnr)
      IF (ALLOCATED (espx)) DEALLOCATE (espx,espy)

      IF (ALLOCATED (kwn)) DEALLOCATE (kwn)
      IF (ALLOCATED (kwnwi)) DEALLOCATE (kwnwi)

      IF (ALLOCATED (elbg)) DEALLOCATE (elbg,relbg,kg)
      IF (ALLOCATED (enr)) DEALLOCATE (enr)
      IF (ALLOCATED (mnr)) DEALLOCATE (mnr)
      IF (ALLOCATED (strnr)) DEALLOCATE (strnr,strom,volt,sigmaa,
     1     kfak,wmatdr,wmatdp,vnr,dat,wmatd,wmatd2,sgmaa2,wdfak)
      IF (ALLOCATED (par)) DEALLOCATE (par,dpar,dpar2)
      IF (ALLOCATED (cgres)) DEALLOCATE (cgres,cgres2)
      IF (ALLOCATED (sigma)) DEALLOCATE (sigma,sigma2)

      IF (ALLOCATED (d0)) DEALLOCATE (d0,fm0)
      IF (ALLOCATED (m0)) DEALLOCATE (m0)

      IF (ALLOCATED (rwddc)) DEALLOCATE (rwddc) 
      IF (ALLOCATED (rwndc)) DEALLOCATE (rwndc) 
      IF (ALLOCATED (rwd)) DEALLOCATE (rwd) 
      IF (ALLOCATED (rwn)) DEALLOCATE (rwn) 
      IF (ALLOCATED (rwdnr)) DEALLOCATE (rwdnr) 

c     Ggf. weiteren Datensatz invertieren
      if (lagain) goto 5

c     'crtomo.cfg' schliessen
      close(fpcfg)

      stop '0'

c.....................................................................

c     Fehlermeldung
 999  open(fperr,file='error.dat',status='replace')
      errflag = 2
      CALL get_error(ftext,errnr,errflag,fetxt)
      write(fperr,'(a80,i3,i1)') fetxt,errnr,errflag
      write(fperr,*)ftext
      close(fperr)
      STOP '-1'

      end

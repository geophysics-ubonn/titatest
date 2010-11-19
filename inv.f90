program inv

!!!$   Hauptprogramm zur Complex-Resistivity-2.5D-Inversion.

!!!$   Belegte Kanaele:  
!!!$   9 - error.dat   -> fperr
!!!$   10 - run.ctr    -> fprun
!!!$   11 - in-/output -> kanal
!!!$   12 - crtomo.cfg -> fpcfg
!!!$   13 - inv.ctr    -> fpinv
!!!$   14 - cjg.ctr    -> fpcjg
!!!$   15 - eps.ctr    -> fpeps
!!!$
!!!$   Andreas Kemna                                        02-May-1995
!!!$   Letzte Aenderung                                     Jul-2010

!!!$.....................................................................

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
  USE invhpmod

!!!$   USE portlib

  IMPLICIT none

  CHARACTER(256)         :: ftext
  INTEGER                :: c1,i
  REAL(KIND(0D0))        :: lamalt
  LOGICAL                :: ols,converged

  INTEGER :: OMP_GET_MAX_THREADS
  INTEGER :: OMP_GET_NUM_THREADS
  INTEGER :: OMP_GET_THREAD_NUM
!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$   SETUP UND INPUT
!!!$   'crtomo.cfg' oeffnen
  fetxt = 'crtomo.cfg'
  errnr = 1
!!!$   Kanal nummern belegen.. die variablen sind global!
  fperr = 9 
  fprun = 10
  kanal = 11 
  fpcfg = 12 
  fpinv = 13 
  fpcjg = 14 
  fpeps = 15 

  open(fpcfg,file=fetxt,status='old',err=999)

  lagain=.TRUE. ! is set afterwards by user input file to false

!!!$   Allgemeine Parameter setzen
  DO WHILE(lagain)
     errnr2 = 0
!!!$   Benoetigte Variablen einlesen
     call rall(kanal,delem,delectr,dstrom,drandb,&
!!!$   diff-     1            dsigma,dvolt,dsens,dstart,lsens,lagain)
!!!$   diff+<
          dsigma,dvolt,dsens,dstart,dd0,dm0,dfm0,lagain)
!!!$   diff+>
     if (errnr.ne.0) goto 999


!!!$   Element- und Randelementbeitraege sowie ggf. Konfigurationsfaktoren
!!!$   zur Berechnung der gemischten Randbedingung bestimmen
     call precal()

     if (errnr.ne.0) goto 999

     if (.not.lbeta) then
        lsr = .false.

!!!$   Ggf. Fehlermeldungen
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

!!!$   Ggf. Fehlermeldungen
        if (swrtr.eq.0) then
           fetxt = ' '
           errnr = 106
        end if
     end if

!!!$   Startmodell belegen
     ALLOCATE (sigma(elanz),sigma2(elanz),stat=errnr)
     IF (errnr /= 0) THEN
        fetxt = 'Error memory allocation fem sigma'
        errnr = 94
        goto 999
     END IF
     call bsigm0(kanal,dstart)
     if (errnr.ne.0) goto 999

!!!$   get memory for model parameters
     ALLOCATE (par(manz),dpar(manz),dpar2(manz),stat=errnr)
     IF (errnr /= 0) THEN
        fetxt = 'Error memory allocation model/update data'
        errnr = 94
        goto 999
     END IF
!!!$   get memory for CG data storage of residuums
     ALLOCATE (cgres(ncgmax+2),cgres2(ncgmax+2),stat=errnr)
     IF (errnr /= 0) THEN
        fetxt = 'Error memory allocation cgres data'
        errnr = 94
        goto 999
     END IF
!!!$   get memory for CJG and update
     ALLOCATE (bvec(manz),STAT=errnr)
     IF (errnr /= 0) THEN
        fetxt = 'Error memory allocation bvec data'
        errnr = 94
        goto 999
     END IF

!!!$   Startparameter setzen
     it     = 0;itr    = 0
     rmsalt = 0d0; lamalt = 1d0; bdpar = 1d0
     betrms = 0d0; pharms = 0d0
     lsetup = .true.; lsetip = .false.; lip    = .false.
     llam   = .false.; ldlami = .true.; lstep  = .false.
     lfstep = .false.
     step   = 1d0; stpalt = 1d0; alam   = 0d0

     CALL wout(kanal,dsigma,dvolt)
!!!$   Kontrolldateien oeffnen
     errnr = 1

     fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
     open(fpinv,file=fetxt,status='replace',err=999)
     close(fpinv)
     fetxt = ramd(1:lnramd)//slash(1:1)//'run.ctr'
     open(fprun,file=fetxt,status='replace',err=999)
!!!$  close(fprun) muss geoeffnet bleiben da sie staendig beschrieben wird
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
        WRITE (fpeps,'(3(G12.5,2x))')&
             (sqrt(wmatd(i)),sqrt(wmatdr(i)),sqrt(wmatdp(i)),i=1,nanz)
     END IF
     close(fpeps)
     errnr = 4

!!!$   Kontrolldateien initialisieren
!!!$   diff-        call kont1(delem,delectr,dstrom,drandb)
!!!$   diff+<
     call kont1(delem,delectr,dstrom,drandb,dd0,dm0,dfm0,lagain)
!!!$   diff+>
     if (errnr.ne.0) goto 999
!!$      CALL SYSTEM('sleep 1000')
!!!$   'sens' zuweisen
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

     write(6,"(a, i3)") " OpenMP max threads: ", OMP_GET_MAX_THREADS()
     !$OMP PARALLEL
     write(6,"(2(a,i3))") " OpenMP: N_threads = ",&
          OMP_GET_NUM_THREADS()," thread = ", OMP_GET_THREAD_NUM()
     !$OMP END PARALLEL

!!!$-------------
!!!$   get current time
     CALL tic(c1)
!!!$.................................................
     converged = .FALSE.
!!!$   MODELLING
!!!$   'a', 'hpot' und 'kpot' zuweisen
     DO WHILE (.NOT.converged)

!!!$   Kontrollausgaben
        WRITE (*,'(a60)',ADVANCE='no')ACHAR(13)//''
        write(*,'(a,i3,a,i3,a)',ADVANCE='no')&
             ACHAR(13)//' Iteration ',it,', ',itr,&
             ' : Calculating Potentials'


        write(fprun,'(a,i3,a,i3,a)')&
             ' Iteration ',it,', ',itr,&
             ' : Calculating Potentials'

        if (ldc) then
           ALLOCATE(adc((mb+1)*sanz),hpotdc(sanz,eanz),bdc(sanz),&
                stat=errnr)
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

        if (ldc) then

!!!$   DC CASE
           do k=1,kwnanz
              do l=1,eanz
                 if (lsr.or.lbeta.or.l.eq.1) then

!!!$   Ggf. Potentialwerte fuer homogenen Fall analytisch berechnen
                    if (lsr) call potana(l,k)

!!!$   Kompilation des Gleichungssystems (fuer Einheitsstrom !)
                    call kompadc(l,k)
                    if (errnr.ne.0) goto 999

!!!$   Ggf. Randbedingung beruecksichtigen
                    if (lrandb) call randdc()
                    if (lrandb2) call randb2()

!!!$   Gleichungssystem skalieren
                    call scaldc()
                    if (errnr.ne.0) goto 999

!!!$   Cholesky-Zerlegung der Matrix
                    call choldc()
                    if (errnr.ne.0) goto 999
                 else

!!!$   Stromvektor modifizieren
                    call kompbdc(l)
                 end if

!!!$   Gleichungssystem loesen
                 call vredc()

!!!$   Potentialwerte zurueckskalieren und umspeichern sowie ggf.
!!!$   analytische Loesung addieren
                 do j=1,sanz
                    kpotdc(j,l,k) = dble(pot(j)) * fak(j)
                    if (lsr) kpotdc(j,l,k) = kpotdc(j,l,k) + &
                         dble(pota(j))
                    if (swrtr.eq.0) hpotdc(j,l) = kpotdc(j,l,k)
                 end do
              end do
           end do

        else

!!!$   COMPLEX CASE
           do k=1,kwnanz
              do l=1,eanz
                 if (lsr.or.lbeta.or.l.eq.1) then

!!!$   Ggf. Potentialwerte fuer homogenen Fall analytisch berechnen
                    if (lsr) call potana(l,k)

!!!$   Kompilation des Gleichungssystems (fuer Einheitsstrom !)
                    call kompab(l,k)
                    if (errnr.ne.0) goto 999

!!!$   Ggf. Randbedingung beruecksichtigen
                    if (lrandb) call randb()
                    if (lrandb2) call randb2()

!!!$   Gleichungssystem skalieren
                    call scalab()
                    if (errnr.ne.0) goto 999

!!!$   Cholesky-Zerlegung der Matrix
                    call chol()
                    if (errnr.ne.0) goto 999
                 else

!!!$   Stromvektor modifizieren
                    call kompb(l)
                 end if

!!!$   Gleichungssystem loesen
                 call vre()

!!!$   Potentialwerte zurueckskalieren und umspeichern sowie ggf.
!!!$   analytische Loesung addieren
                 do j=1,sanz
                    kpot(j,l,k) = pot(j) * dcmplx(fak(j))
                    if (lsr) kpot(j,l,k) = kpot(j,l,k) + pota(j)
                    if (swrtr.eq.0) hpot(j,l) = kpot(j,l,k)
                 end do
              end do
           end do

        end if

!!!$   Ggf. Ruecktransformation der Potentialwerte
        if (swrtr.eq.1) call rtrafo()

!!!$   Spannungswerte berechnen
        call bvolti()
        if (errnr.ne.0) goto 999

!!!$   'a' und 'hpot' freigeben
        if (ldc) then
           DEALLOCATE(adc,hpotdc,bdc)
        else
           DEALLOCATE(a,hpot,b)
        end if
        if (lsetup) then

!!!$   Ggf. background auf ratio-Daten "multiplizieren"
           if (lratio) then
              do j=1,nanz
                 dat(j) = dat(j) + sigmaa(j)
              end do
           end if

!!!$   Polaritaeten checken
           call chkpol(lsetup)
        end if

!!!$   Daten-RMS berechnen
        call dmisft(lsetup)

        if (errnr.ne.0) goto 999

!!!$   'nrmsd=0' ausschliessen
        if (nrmsd.lt.1d-12) nrmsd=nrmsdm*(1d0-mqrms)

!!!$   tst
!!!$   tst        if (lfphai) then
!!!$   tst            llam = .true.
!!!$   tst            if (.not.lip) nrmsd = 1d0
!!!$   tst        end if

!!!$.............................

!!!$   Kontrollvariablen ausgeben
        call kont2(lsetup)
        if (errnr.ne.0) goto 999

!!!$   ABBRUCHBEDINGUNGEN
        if (llam.and..not.lstep) then
           lamalt = lam 
!!!$   Polaritaeten checken
           call chkpol(lsetup)

!!!$   Wiederholt minimale step-length ?
           if (stpalt.eq.0d0) errnr2=92

!!!$   Keine Verbesserung des Daten-RMS ?
           if (dabs(1d0-rmsalt/nrmsd).le.mqrms) errnr2=81

!!!$   Minimaler Daten-RMS erreicht ?
!!!$   tst            if (dabs(1d0-nrmsd/nrmsdm).le.mqrms) errnr2=80
           if (dabs(1d0-nrmsd/nrmsdm).le.mqrms.and.ldlamf) errnr2=80

!!!$   Maximale Anzahl an Iterationen ?
           if (it.ge.itmax) errnr2=79

!!!$   Minimal stepsize erreicht ?
           IF (bdpar < bdmin) THEN
              errnr2=109
              WRITE (ftext,*)'check stepsize',bdpar,it,itr
           END IF

!!!$   Ggf. abbrechen oder "final phase improvement"
           if (errnr2.ne.0) then
              if (lfphai.and.errnr2.ne.79) then
                 errnr2 = 0
!!!$   ak
!!!$   Widerstandsverteilung und modellierte Daten ausgeben
                 call wout(kanal,dsigma,dvolt)
                 if (errnr.ne.0) goto 999

!!!$   Kontrollausgaben
                 write(*,'(/a/)')&
                      '**** Final phase improvement ****'

                 write(fprun,'(a24)',err=999)&
                      ' Final phase improvement'

                 fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
                 open(fpinv,file=fetxt,status='old',POSITION='append',err=999)
                 write(fpinv,'(/a/)',err=999)&
                      '------------------------------------------------'//&
                      '------------------------------------------------'//&
                      '-----------------'
                 close(fpinv)

!!!$   Wichtungsfeld umspeichern
                 wmatd = wmatdp

                 lip    = .true.
                 lsetip = .true. ! 
                 lfphai = .false.
                 llam   = .false.
                 ldlami = .true.
                 lfstep = .true.
                 step   = 1d0

!!!$   ak
                 fetxt = 'cp -f inv.lastmod inv.lastmod_rho'
                 CALL SYSTEM (TRIM(fetxt))
                 IF (lffhom) THEN
                    write(*,*)&
                         ' ******* Restarting phase model ********'
                    write(fprun,*)&
                         ' ******* Restarting phase model ********'
                    do j=1,elanz
                       sigma(j) = dcmplx(&
                            dcos(pha0/1d3)*cdabs(sigma(j)) ,&
                            -dsin(pha0/1d3)*cdabs(sigma(j)) )
                    end do

                    CYCLE       ! neues calc
                 END IF
!!!$   ak

!!!$   Daten-RMS berechnen
                 call dmisft(lsetip)
                 if (errnr.ne.0) goto 999

!!!$   Kontrollvariablen ausgeben
                 call kont2(lsetip)
                 if (errnr.ne.0) goto 999
              else
                 EXIT
              end if
           else
!!!$   ak
!!!$   Widerstandsverteilung und modellierte Daten ausgeben
              call wout(kanal,dsigma,dvolt)
              if (errnr.ne.0) goto 999
           end if
        end if

        IF (ldc .OR. lip) THEN
           ALLOCATE (pvecdc(manz),bvecdc(manz),rvecdc(manz), &
                apdc(nanz),STAT=errnr)
           IF (errnr /= 0) THEN
              fetxt = 'Error memory allocation rvec/bvec dc'
              errnr = 94
              GOTO 999
           END IF
        ELSE
           ALLOCATE (rvec(manz),pvec(manz),ap(nanz),stat=errnr)
           IF (errnr /= 0) THEN
              fetxt = 'Error memory allocation rve!!!$in cjggdc'
              errnr = 94
              GOTO 999
           END IF
        END IF


        if ((llam.and..not.lstep).or.lsetup.or.lsetip) then
!!!$   Iterationsindex hochzaehlen
           it = it+1

!!!$   Parameter zuruecksetzen
           itr    = 0
           rmsreg = 0d0
           dlam   = 1d0
           dlalt  = 1d0
           ldlamf = .false.

!!!$   Daten-RMS speichern
           rmsalt = nrmsd

!!!$   Lambda speichen
           IF (it>1) lamalt = lam ! mindestens einmal konvergiert

!!!$   Felder speichern
           sigma2 = sigma

           sgmaa2 = sigmaa

           IF (lrobust) wmatd2 = wmatd

!!!$   Kontrollausgaben
           write (*,'(a,i3,a,i3,a)',ADVANCE='no') &
                ACHAR(13)//' Iteration ',it,', ',itr,&
                ' : Calculating Sensitivities'

           write(fprun,'(a,i3,a,i3,a)')' Iteration ',it,', ',itr,&
                ' : Calculating Sensitivities'

!!!$   SENSITIVITAETEN berechnen
           if (ldc) then
              call bsendc()
           else
              call bsensi()
           end if

!!!$   evtl   Rauhigkeitsmatrix belegen
           IF (lsetup.OR.(ltri>3.AND.ltri<10)) CALL bsmatm(it) 

        else
!!!$   Felder zuruecksetzen
           sigma = sigma2

           sigmaa = sgmaa2

           IF (lrobust) wmatd = wmatd2

        end if
        IF (itmax == 0) THEN ! only precalcs..
           errnr2 = 109
           print*,'Only precalcs'
           EXIT
        END IF
!!!$!!!$   'kpot' freigeben wird nicht mehr freigegeben
!!!$      if (ldc) then
!!!$         DEALLOCATE(kpotdc)
!!!$      else
!!!$         DEALLOCATE(kpot)
!!!$      end if

!!!$   REGULARISIERUNG / STEP-LENGTH einstellen
        if (.not.lstep) then
           if (llam) then

!!!$   "Regularisierungsschleife" initialisieren und step-length zuruecksetzen
              llam = .false.
              step = 1d0
           else

!!!$   Regularisierungsindex hochzaehlen
              itr = itr+1
              if (((((nrmsd.lt.rmsreg.and.itr.le.nlam).or. &
                   (dlam.gt.1d0.and.itr.le.nlam)).and.&
                   (.not.ldlamf.or.dlalt.le.1d0).and.&
                   dabs(1d0-rmsreg/nrmsdm).gt.mqrms).or.&
                   (rmsreg.eq.0d0)).AND.&
                   (bdpar >= bdmin)) then
!!!$   Regularisierungsparameter bestimmen
                 if (lsetup.or.lsetip) then

!!!$   Kontrollausgabe
                    IF (llamf) THEN
                       lam = lamfix
                    ELSE
                       write(*,'(a,i3,a,i3,a)',ADVANCE='no')&
                            ACHAR(13)//' Iteration ',it,', ',itr,&
                            ' : Calculating 1st regularization parameter'
                       write(fprun,'(a,i3,a,i3,a)',ADVANCE='no')&
                            ' Iteration ',it,', ',itr,&
                            ' : Calculating 1st regularization parameter'
                       call blam0()
                       WRITE (*,'(a,G10.2)',ADVANCE='no')'lam_0:: ',lammax
                       WRITE (fprun,'(a,G10.2)')'lam_0 ',lammax
                       lam = lammax
!!!$   ak Model EGS2003, ERT2003                        call blam0()
!!!$   ak Model EGS2003, ERT2003                        lam = lammax
!!!$   ak                        lam = 1d4
                    END IF
                 else
                    IF (llamf) THEN
                       lam = lamfix
                    ELSE
                       dlalt = dlam
                       if (ldlami) then
                          ldlami = .false.
                          alam   = dmax1(dabs(dlog(nrmsd/nrmsdm)),&
                               dlog(1d0+mqrms))
                          dlam   = fstart
                       else
                          alam = dmax1(alam,dabs(dlog(nrmsd/nrmsdm)))
                          dlam = dlog(fstop)*&
                               sign(1d0,dlog(nrmsd/nrmsdm))+&
                               dlog(fstart/fstop)*&
                               dlog(nrmsd/nrmsdm)/alam
                          dlam = dexp(dlam)
                       end if
                       lam = lam*dlam
                       if (dlalt.gt.1d0.and.dlam.lt.1d0) ldlamf=.true.
!!!$   tst                        if (dlam.gt.1d0) lfstep=.true.
!!!$   ak Model EGS2003
                       if (dlam.gt.1d0) lrobust=.false.
                    END IF
                 end if
              else

!!!$   Regularisierungsparameter zuruecksetzen und step-length verkleinern
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

!!!$   Ggf. Daten-RMS speichern
              if (lsetup.or.lsetip) then
                 lsetup = .false.
                 lsetip = .false.
              else
                 if (.not.lstep) rmsreg=nrmsd
              end if
           end if
        else
           lstep = .false.

!!!$   Parabolische Interpolation zur Bestimmung der optimalen step-length
           call parfit(rmsalt,nrmsd,rmsreg,nrmsdm,stpmin)

           if (step.eq.stpmin.and.stpalt.eq.stpmin)then

!!!$   Nach naechstem Modelling abbrechen
              stpalt = 0d0
           else

!!!$   Step-length speichern
              stpalt = step
           end if
        end if
!!!$   Kontrollausgaben
        write(*,'(a,i3,a,i3,a)',ADVANCE='no')&
             ACHAR(13)//' Iteration ',it,', ',itr,&
             ' : Updating'

        write(fprun,*)' Iteration ',it,', ',itr,&
             ' : Updating'

!!!$ Modell parameter mit aktuellen Leitfaehigkeiten belegen

        CALL bpar
        if (errnr.ne.0) goto 999

!!!$   UPDATE anbringen
        call update
        if (errnr.ne.0) goto 999

        IF (ALLOCATED (rvecdc)) &
             DEALLOCATE (rvecdc,apdc,bvecdc,pvecdc,STAT=errnr)

        IF (ALLOCATED (pvec)) &
             DEALLOCATE (rvec,pvec,ap,STAT=errnr)

!!!$ Leitfaehigkeiten mit verbessertem Modell belegen
        CALL bsigma
        if (errnr.ne.0) goto 999

!!!$   Roughness bestimmen
        CALL brough

!!!$   Ggf. Referenzleitfaehigkeit bestimmen
        if (lsr) call refsig()
!!!$   Neues Modelling
     END DO ! while not converged

!!!$.................................................

!!!$   OUTPUT
     call wout(kanal,dsigma,dvolt)
     if (errnr.ne.0) goto 999

!!!$   Ggf. Summe der Sensitivitaeten aller Messungen ausgeben
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

!!!$   Kontrollausgaben

     if (errnr2.eq.92) then

        write(*,'(a22,a31)') ' Iteration terminated:',&
             ' Min. step-length for 2nd time.'

        write(fprun,'(a22,a31)',err=999) ' Iteration terminated:',&
             ' Min. step-length for 2nd time.'
     else if (errnr2.eq.80) then
        write(*,'(a22,a10)') ' Iteration terminated:', &
             ' Min. RMS.'

        write(fprun,'(a22,a10)',err=999) ' Iteration terminated:',&
             ' Min. RMS.'
     else if (errnr2.eq.81) then
        write(*,'(a22,a24)') ' Iteration terminated:',&
             ' Min. rel. RMS decrease.'

        write(fprun,'(a22,a24)',err=999) ' Iteration terminated:',&
             ' Min. rel. RMS decrease.'
     else if (errnr2.eq.79) then
        write(*,'(a22,a19)') ' Iteration terminated:',&
             ' Max. # iterations.'

        write(fprun,'(a22,a19)',err=999) ' Iteration terminated:',&
             ' Max. # iterations.'

     else if (errnr2.eq.109) then
        write(*,'(a)') ' Iteration terminated:'//&
             ' Min. model changes reached'

        write(fprun,'(a)',err=999) ' Iteration terminated:'//&
             ' Min. model changes reached'

     end if

!!!$   Run-time abfragen und ausgeben
     fetxt = ' CPU time: '
     CALL toc(c1,fetxt)
     write(fprun,'(a)',err=999)TRIM(fetxt)

!!!$   Kontrolldateien schliessen
     close(fprun)
     close(fpinv)
     close(fpcjg)
     close(fpeps)

!!!$   'sens' und 'kpot' freigeben
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
     IF (ALLOCATED (strnr)) DEALLOCATE (strnr,strom,volt,sigmaa,&
          kfak,wmatdr,wmatdp,vnr,dat,wmatd,wmatd2,sgmaa2,wdfak)
     IF (ALLOCATED (par)) DEALLOCATE (par,dpar,dpar2)
     IF (ALLOCATED (cgres)) DEALLOCATE (cgres,cgres2)

     IF (ALLOCATED (bvec)) DEALLOCATE (bvec)

     IF (ALLOCATED (sigma)) DEALLOCATE (sigma,sigma2)

     IF (ALLOCATED (d0)) DEALLOCATE (d0,fm0)
     IF (ALLOCATED (m0)) DEALLOCATE (m0)

     IF (ALLOCATED (rwddc)) DEALLOCATE (rwddc) 
     IF (ALLOCATED (rwndc)) DEALLOCATE (rwndc) 
     IF (ALLOCATED (rwd)) DEALLOCATE (rwd) 
     IF (ALLOCATED (rwn)) DEALLOCATE (rwn) 
     IF (ALLOCATED (rwdnr)) DEALLOCATE (rwdnr) 

!!!$   Ggf. weiteren Datensatz invertieren
  END DO

!!!$   'crtomo.cfg' schliessen
  close(fpcfg)

  stop '0'

!!!$.....................................................................

!!!$   Fehlermeldung
999 open(fperr,file='error.dat',status='replace')
  errflag = 2
  CALL get_error(ftext,errnr,errflag,fetxt)
  write(fperr,'(a80,i3,i1)') fetxt,errnr,errflag
  write(fperr,*)ftext
  close(fperr)
  STOP '-1'

end program inv

PROGRAM inv

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
  USE omp_lib
  USE ompmod
  USE get_ver
!!!$   USE portlib

  IMPLICIT NONE

  CHARACTER(256)         :: ftext
  INTEGER                :: c1,i,count,mythreads,maxthreads
  REAL(KIND(0D0))        :: lamalt
  LOGICAL                :: converged,l_bsmat

  INTEGER :: getpid,pid
!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$   SETUP UND INPUT
!!!$   'crtomo.cfg' oeffnen
  errnr = 1
!!!$   Kanal nummern belegen.. die variablen sind global!
  fperr = 9 
  fprun = 10
  kanal = 11 
  fpcfg = 12 
  fpinv = 13 
  fpcjg = 14 
  fpeps = 15 

  pid = getpid()
  fetxt = 'crtomo.pid'
  PRINT*,'CRTomo Process_ID ::',pid
  OPEN (fprun,FILE=TRIM(fetxt),STATUS='replace',err=999)
  WRITE (fprun,*)pid
  CLOSE (fprun)
  maxthreads = OMP_GET_MAX_THREADS()
  WRITE(6,"(a, i3)") " OpenMP max threads: ", maxthreads

  CALL get_git_ver

  fetxt = 'crtomo.cfg'
  OPEN(fpcfg,file=TRIM(fetxt),status='old',err=999)

  lagain=.TRUE. ! is set afterwards by user input file to false

!!!$   Allgemeine Parameter setzen
  DO WHILE ( lagain ) ! this loop exits after all files are processed

     errnr2 = 0
!!!$   Benoetigte Variablen einlesen
     CALL rall(kanal,delem,delectr,dstrom,drandb,&
          dsigma,dvolt,dsens,dstart,dd0,dm0,dfm0,lagain)
!!!$   diff+<
!!!$   diff-     1            dsigma,dvolt,dsens,dstart,lsens,lagain)
!!!$   diff+>
     IF (errnr.NE.0) GOTO 999

     NTHREADS = maxthreads
!!!$ now that we know nf and kwnanz, we can adjust the OMP environment..
     IF (maxthreads > 2) THEN ! single or double processor machines don't need scheduling..
        mythreads = MAX(kwnanz,2)
        PRINT*,'Rescheduling..'
        IF ( mythreads <= maxthreads ) THEN ! best case,
!!!$ the number of processors is greater or equal the assumed
!!!$ workload
           PRINT*,'perfect match'
        ELSE 
!!!$ is smaller than the minimum workload.. now we have to devide a bit..
           PRINT*,'less nodes than wavenumbers'
           DO i = 1, INT(kwnanz/2)
              mythreads = INT(kwnanz / i) + 1
              IF (mythreads < maxthreads) EXIT
           END DO
        END IF
        NTHREADS = mythreads
     END IF
     CALL OMP_SET_NUM_THREADS ( NTHREADS )
     ! recheck ..
     i = OMP_GET_MAX_THREADS()
     WRITE(6,'(2(a, i3),a)') " OpenMP threads: ",i,'(',maxthreads,')'
!!!$
     
!!!$   Element- und Randelementbeitraege sowie ggf. Konfigurationsfaktoren
!!!$   zur Berechnung der gemischten Randbedingung bestimmen
     CALL precal()

     IF (errnr.NE.0) GOTO 999

     IF (.NOT.lbeta) THEN
        lsr = .FALSE.

!!!$   Ggf. Fehlermeldungen
        IF (.NOT.lsink) THEN
           fetxt = 'no mixed boundary specify sink node'
           errnr = 102
        END IF
!!!$ RM this is handeled in rall..
!!$        if (swrtr.eq.0.and..not.lrandb2) then
!!$           fetxt = ' '
!!$           errnr = 103
!!$        end if

        IF (errnr.NE.0) GOTO 999
     ELSE

!!!$   Ggf. Fehlermeldungen
        IF (swrtr.EQ.0) THEN
           fetxt = ' '
           errnr = 106
        END IF
     END IF


!!!$   getting dynamic memory 
     errnr = 94
!!!$ physical model
     fetxt = 'allocation problem sigma'
     ALLOCATE (sigma(elanz),STAT=errnr)
     IF (errnr /= 0) GOTO 999
     fetxt = 'allocation problem sigma2'
     ALLOCATE (sigma2(elanz),STAT=errnr)
     IF (errnr /= 0) GOTO 999
!!!$  model parameters
     fetxt = 'allocation problem par'
     ALLOCATE (par(manz),STAT=errnr)
     IF (errnr /= 0) GOTO 999
     fetxt = 'allocation problem dpar'
     ALLOCATE (dpar(manz),STAT=errnr)
     IF (errnr /= 0) GOTO 999
     fetxt = 'allocation problem dpar2'
     ALLOCATE (dpar2(manz),STAT=errnr)
     IF (errnr /= 0) GOTO 999
     fetxt = 'allocation problem pot'
     ALLOCATE(pot(sanz),STAT=errnr)
     IF (errnr /= 0) GOTO 999
     fetxt = 'allocation problem pota'
     ALLOCATE (pota(sanz),STAT=errnr)
     IF (errnr /= 0) GOTO 999
!!!$ now the big array are coming.. 
     fetxt = 'allocation problem fak'
     ALLOCATE (fak(sanz),STAT=errnr) ! fak for modeling
     IF (errnr /= 0) GOTO 999
     IF (ldc) THEN
        fetxt = 'allocation problem sensdc'
        ALLOCATE (sensdc(nanz,manz),STAT=errnr)
        IF (errnr /= 0) GOTO 999
        fetxt = 'allocation problem kpotdc'
        ALLOCATE (kpotdc(sanz,eanz,kwnanz),STAT=errnr)
     ELSE
        fetxt = 'allocation problem sens'
        ALLOCATE (sens(nanz,manz),STAT=errnr)
        IF (errnr /= 0) GOTO 999
        fetxt = 'allocation problem kpot'
        ALLOCATE (kpot(sanz,eanz,kwnanz),STAT=errnr)
     END IF
     IF (errnr /= 0) GOTO 999


!!!$ get CG data storage of residuums and bvec, which is global
     CALL con_cjgmod (1,fetxt,errnr)
     IF (errnr /= 0) GOTO 999


!!!$ set starting model 
     CALL bsigm0(kanal,dstart)
     IF (errnr.NE.0) GOTO 999

!!!$   Startparameter setzen
     it     = 0;itr    = 0
     rmsalt = 0d0; lamalt = 1d0; bdpar = 1d0
     IF (llamf) lamalt = lamfix
     betrms = 0d0; pharms = 0d0
     lsetup = .TRUE.; lsetip = .FALSE.; lip    = .FALSE.
     llam   = .FALSE.; ldlami = .TRUE.; lstep  = .FALSE.
     lfstep = .FALSE.; l_bsmat = .TRUE.
     step   = 1d0; stpalt = 1d0; alam   = 0d0

!!!$   Kontrolldateien oeffnen
     errnr = 1

     fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
     OPEN(fpinv,file=TRIM(fetxt),status='replace',err=999)
     CLOSE(fpinv)
     fetxt = ramd(1:lnramd)//slash(1:1)//'run.ctr'
     OPEN(fprun,file=TRIM(fetxt),status='replace',err=999)
!!!$  close(fprun) muss geoeffnet bleiben da sie staendig beschrieben wird
     fetxt = ramd(1:lnramd)//slash(1:1)//'cjg.ctr'
     OPEN(fpcjg,file=TRIM(fetxt),status='replace',err=999)
     CLOSE(fpcjg)
     fetxt = ramd(1:lnramd)//slash(1:1)//'eps.ctr'
     OPEN(fpeps,file=TRIM(fetxt),status='replace',err=999)

!!!$  Write errors for all measurements to fpeps
     IF (ldc) THEN
        WRITE (fpeps,'(a)')'1/eps_r      datum'
        WRITE (fpeps,'(G12.3,2x,G14.5)')(SQRT(wmatdr(i)),&
             REAL(dat(i)),i=1,nanz)
     ELSE
        WRITE (fpeps,'(t5,a,t14,a,t27,a,t38,a,t50,a,t62,a,t71,a,t87,a)')'1/eps_r','1/eps_p',&
             '1/eps','eps_r','eps_p','eps','-log(|R|)', '-Phase (rad)'
        WRITE (fpeps,'(3F10.1,2x,3G12.3,2G15.7)')&
             (SQRT(wmatdr(i)),SQRT(wmatdp(i)),SQRT(wmatd(i)),1/SQRT(wmatdr(i)),&
             1/SQRT(wmatdp(i)),1/SQRT(wmatd(i)),REAL(dat(i)),AIMAG(dat(i)),i=1,nanz)
     END IF
     CLOSE(fpeps)
     errnr = 4

!!!$   Kontrolldateien initialisieren
!!!$   diff-        call kont1(delem,delectr,dstrom,drandb)
!!!$   diff+<
     CALL kont1(delem,delectr,dstrom,drandb,dd0,dm0,dfm0,lagain)
!!!$   diff+>
     IF (errnr.NE.0) GOTO 999

!!$     !$OMP PARALLEL
!!$     write(6,"(2(a,i3))") " OpenMP: N_threads = ",&
!!$          OMP_GET_NUM_THREADS()," thread = ", OMP_GET_THREAD_NUM()
!!$     !$OMP END PARALLEL

!!!$-------------
!!!$   get current time
     CALL tic(c1)
!!!$.................................................
     converged = .FALSE.

     DO WHILE (.NOT.converged) ! optimization loop

!!!$   Control output
        WRITE(*,'(a,i3,a,i3,a,t100,a)',ADVANCE='no')ACHAR(13)//&
             ' Iteration ',it,', ',itr,' : Calculating Potentials',''
        WRITE(fprun,'(a,i3,a,i3,a)')' Iteration ',it,', ',itr,&
             ' : Calculating Potentials'

!!!$   MODELLING
        count = 0
        IF (ldc) THEN
           fetxt = 'allocation problem adc'
           ALLOCATE (adc((mb+1)*sanz),STAT=errnr)
           IF (errnr /= 0) GOTO 999
           fetxt = 'allocation problem hpotdc'
           ALLOCATE (hpotdc(sanz,eanz),STAT=errnr)
           IF (errnr /= 0) GOTO 999
           ALLOCATE (bdc(sanz),STAT=errnr)
           fetxt = 'allocation problem adc'
           IF (errnr /= 0) GOTO 999

!$OMP PARALLEL DEFAULT (none) &
!$OMP FIRSTPRIVATE (pota,fak,pot,adc,bdc,fetxt) &
!$OMP PRIVATE (j,l) &
!$OMP SHARED (kwnanz,lverb,eanz,lsr,lbeta,lrandb,lrandb2,sanz,kpotdc,swrtr,hpotdc,elbg,count)           
!$OMP DO 
!!!$   DC CASE
           DO k=1,kwnanz
              !$OMP ATOMIC
              count = count + 1
              fetxt = 'DC-Calculation wavenumber'
              IF (lverb) WRITE (*,'(a,t35,I4,t100,a)',ADVANCE='no')&
                   ACHAR(13)//TRIM(fetxt),count,''
              DO l=1,eanz
                 IF (lsr.OR.lbeta.OR.l.EQ.1) THEN
!!!$   Evtl calculation of analytical potentials
                    IF (lsr) CALL potana(l,k,pota)

!!!$   COMPilation of the linear system
                    fetxt = 'kompadc'
                    CALL kompadc(l,k,adc,bdc)
!                    if (errnr.ne.0) goto 999

!!!$   Evtl take Dirichlet boundary values into account
                    IF (lrandb) CALL randdc(adc,bdc)
                    IF (lrandb2) CALL randbdc2(adc,bdc)

!!!$   Scale the linear system (preconditioning stores fak)
                    fetxt = 'scaldc'
                    CALL scaldc(adc,bdc,fak)
!                    if (errnr.ne.0) goto 999
!!!$   Cholesky-Factorization of the Matrix
                    fetxt = 'choldc'
                    CALL choldc(adc)
!                    if (errnr.ne.0) goto 999

                 ELSE
                    fetxt = 'kompbdc'
!!!$   Modification of the current vector (Right Hand Side)
                    CALL kompbdc(l,bdc,fak)
                 END IF

!!!$   Solve linear system
                 fetxt = 'vredc'
                 CALL vredc(adc,bdc,pot)
!!!$   Scale back the potentials, save them and 
!!!$   eventually add the analytical response
                 DO j=1,sanz
                    kpotdc(j,l,k) = DBLE(pot(j)) * fak(j)
                    IF (lsr) kpotdc(j,l,k) = kpotdc(j,l,k) + &
                         DBLE(pota(j))
                    IF (swrtr.EQ.0) hpotdc(j,l) = kpotdc(j,l,k)
                 END DO
              END DO
           END DO
!$OMP END DO
!$OMP END PARALLEL

        ELSE

           fetxt = 'allocation problem a'
           ALLOCATE (a((mb+1)*sanz),STAT=errnr)
           IF (errnr /= 0) GOTO 999
           fetxt = 'allocation problem hpot'
           ALLOCATE (hpot(sanz,eanz),STAT=errnr) 
           IF (errnr /= 0) GOTO 999
           fetxt = 'allocation problem b'
           ALLOCATE (b(sanz),STAT=errnr)
           IF (errnr /= 0) GOTO 999
!$OMP PARALLEL DEFAULT (none) &
!$OMP FIRSTPRIVATE (pota,fak,pot,a,b,fetxt) &
!$OMP PRIVATE (j,l,k) &
!$OMP SHARED (kwnanz,lverb,eanz,lsr,lbeta,lrandb,lrandb2,sanz,kpot,swrtr,hpot,count)
           !$OMP DO
!!!$   COMPLEX CASE
           DO k=1,kwnanz
              !$OMP ATOMIC
              count = count + 1
              fetxt = 'IP-Calculation wavenumber'
              IF (lverb) WRITE (*,'(a,t35,I4,t100,a)',ADVANCE='no')&
                   ACHAR(13)//TRIM(fetxt),count,''
              DO l=1,eanz
                 IF (lsr.OR.lbeta.OR.l.EQ.1) THEN

!!!$   Ggf. Potentialwerte fuer homogenen Fall analytisch berechnen
                    IF (lsr) CALL potana(l,k,pota)

!!!$   KOMPilation des Gleichungssystems (fuer Einheitsstrom !)
                    fetxt = 'kompab'
                    CALL kompab(l,k,a,b)
!                    if (errnr.ne.0) goto 999

!!!$   Ggf. Randbedingung beruecksichtigen
                    IF (lrandb) CALL randb(a,b)
                    IF (lrandb2) CALL randb2(a,b)

!!!$   Gleichungssystem skalieren
                    fetxt = 'scalab'
                    CALL scalab(a,b,fak)
!                    if (errnr.ne.0) goto 999

!!!$   Cholesky-Zerlegung der Matrix
                    fetxt = 'chol'
                    CALL chol(a)
!                    if (errnr.ne.0) goto 999
                 ELSE

!!!$   Stromvektor modifizieren
                    fetxt = 'kompb'
                    CALL kompb(l,b,fak)

                 END IF

!!!$   Gleichungssystem loesen
                 fetxt = 'vre'
                 CALL vre(a,b,pot)

!!!$   Potentialwerte zurueckskalieren und umspeichern sowie ggf.
!!!$   analytische Loesung addieren
                 DO j=1,sanz
                    kpot(j,l,k) = pot(j) * dcmplx(fak(j))
                    IF (lsr) kpot(j,l,k) = kpot(j,l,k) + pota(j)
                    IF (swrtr.EQ.0) hpot(j,l) = kpot(j,l,k)
                 END DO
              END DO
           END DO
!$OMP END DO
!$OMP END PARALLEL

        END IF

!!!$   Ggf. Ruecktransformation der Potentialwerte
        IF (swrtr.EQ.1) THEN
           CALL rtrafo(errnr)
           IF (errnr.NE.0) GOTO 999
        END IF
!!!$   Spannungswerte berechnen
        CALL bvolti()
        IF (errnr.NE.0) GOTO 999
!!$  free some memory..
        IF (ldc) THEN
           DEALLOCATE(adc,hpotdc,bdc)
        ELSE
           DEALLOCATE(a,hpot,b)
        END IF

        IF (lsetup.OR.lsetip) THEN
!!!$   Ggf. background auf ratio-Daten "multiplizieren"
           IF (lratio) THEN
              DO j=1,nanz
                 dat(j) = dat(j) + sigmaa(j)
              END DO
           END IF

!!!$   Polaritaeten checken
           CALL chkpol(lsetup.OR.lsetip)
        END IF

!!!$   Daten-RMS berechnen
        CALL dmisft(lsetup.OR.lsetip)
!        print*,nrmsd,betrms,pharms,lrobust,l1rat
        IF (errnr.NE.0) GOTO 999
        WRITE (*,'(a,F8.3)',ADVANCE='no')'actual fit',nrmsd
!!!$   'nrmsd=0' ausschliessen
        IF (nrmsd.LT.1d-12) nrmsd=nrmsdm*(1d0-mqrms)

!!!$   tst
!!!$   tst        if (lfphai) then
!!!$   tst            llam = .true.
!!!$   tst            if (.not.lip) nrmsd = 1d0
!!!$   tst        end if

!!!$.............................
        IF (it == 0) THEN
           WRITE (*,'(/a,t100/)')ACHAR(13)//&
                'WRITING STARTING MODEL'
           CALL wout(kanal,dsigma,dvolt)
        END IF
!!!$   Kontrollvariablen ausgeben
        CALL kont2(lsetup.OR.lsetip)
        IF (errnr.NE.0) GOTO 999

!!!$   ABBRUCHBEDINGUNGEN
        IF (llam.AND..NOT.lstep) THEN
           lamalt = lam 
!!!$   Polaritaeten checken
           CALL chkpol(lsetup.OR.lsetip)

!!!$   Wiederholt minimale step-length ?
           IF (stpalt.EQ.0d0) errnr2=92

!!!$   Keine Verbesserung des Daten-RMS ?
           IF (dabs(1d0-rmsalt/nrmsd).LE.mqrms) errnr2=81

!!!$   Minimaler Daten-RMS erreicht ?
!!!$   tst            if (dabs(1d0-nrmsd/nrmsdm).le.mqrms) errnr2=80
           IF (dabs(1d0-nrmsd/nrmsdm).LE.mqrms.AND.ldlamf) errnr2=80

!!!$   Maximale Anzahl an Iterationen ?
           IF (it.GE.itmax) errnr2=79

!!!$   Minimal stepsize erreicht ?
           IF (bdpar < bdmin) THEN
              errnr2=109
              WRITE (ftext,*)'check stepsize',bdpar,it,itr
           END IF

!!!$   Ggf. abbrechen oder "final phase improvement"
           IF (errnr2.NE.0) THEN
              IF (lfphai.AND.errnr2.NE.79) THEN
                 errnr2 = 0
!!!$   ak
!!!$   Widerstandsverteilung und modellierte Daten ausgeben
                 CALL wout(kanal,dsigma,dvolt)
                 IF (errnr.NE.0) GOTO 999

!!!$   Kontrollausgaben
                 WRITE(*,'(/a/)')&
                      '**** Final phase improvement ****'

                 WRITE(fprun,'(a24)',err=999)&
                      ' Final phase improvement'

                 fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
                 OPEN(fpinv,file=TRIM(fetxt),status='old',&
                      POSITION='append',err=999)
                 WRITE(fpinv,'(/a/)',err=999)&
                      '------------------------------------------------'//&
                      '------------------------------------------------'//&
                      '-----------------'
                 CLOSE(fpinv)

!!!$   Wichtungsfeld umspeichern
                 wmatd = wmatdp

                 lip    = .TRUE.
                 lsetip = .TRUE. ! 
                 lfphai = .FALSE.
                 llam   = .FALSE.
                 ldlami = .TRUE.
                 lfstep = .TRUE.
                 step   = 1d0

!!!$   ak
                 fetxt = 'cp -f inv.lastmod inv.lastmod_rho'
                 CALL SYSTEM (TRIM(fetxt))
                 IF (lffhom) THEN
                    WRITE(*,*)&
                         ' ******* Restarting phase model ********'
                    WRITE(fprun,*)&
                         ' ******* Restarting phase model ********'
                    fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
                    OPEN (fpinv,file=TRIM(fetxt),status='old',err=999,&
                         position='append')
                    WRITE (fpinv,*)&
                         ' ******* Resetting phase model ********'
                    CLOSE (fpinv)
                    DO j=1,elanz
                       sigma(j) = dcmplx(&
                            dcos(pha0/1d3)*cdabs(sigma(j)) ,&
                            -dsin(pha0/1d3)*cdabs(sigma(j)) )
                    END DO
                    lsetup = .TRUE. ! ensure proper misfit and kont2 output
                    CYCLE       ! neues calc

                 END IF
!!!$   ak

!!!$   Daten-RMS berechnen
                 CALL dmisft(lsetip)
                 IF (errnr.NE.0) GOTO 999

!!!$   Kontrollvariablen ausgeben
                 CALL kont2(lsetip)
                 IF (errnr.NE.0) GOTO 999
              ELSE
                 EXIT
              END IF
           ELSE
!!!$   ak
!!!$   Widerstandsverteilung und modellierte Daten ausgeben
!!$              WRITE (*,'(a,t30,I4,t100,a)')ACHAR(13)//&
!!$                   'WRITING MODEL ITERATE',it,''
              CALL wout(kanal,dsigma,dvolt)
              IF (errnr.NE.0) GOTO 999
           END IF
        END IF

        IF ((llam.AND..NOT.lstep).OR.lsetup.OR.lsetip) THEN

!!!$   Iterationsindex hochzaehlen
           it = it+1

!!!$   Parameter zuruecksetzen
           itr    = 0
           rmsreg = 0d0
           dlam   = 1d0
           dlalt  = 1d0
           ldlamf = .FALSE.

!!!$   Daten-RMS speichern
           rmsalt = nrmsd

!!!$   Lambda speichen
           IF (it>1) lamalt = lam ! mindestens einmal konvergiert

!!!$   Felder speichern
           sigma2 = sigma

           sgmaa2 = sigmaa

           IF (lrobust) wmatd2 = wmatd

!!!$   Kontrollausgaben
           WRITE (*,'(a,i3,a,i3,a,t63,a,t78,F8.3,t100,a)',ADVANCE='no') &
                ACHAR(13)//' Iteration ',it,', ',itr,&
                ' : Calculating Sensitivities','fit',nrmsd,''

           WRITE(fprun,'(a,i3,a,i3,a)')' Iteration ',it,', ',itr,&
                ' : Calculating Sensitivities'

!!!$   SENSITIVITAETEN berechnen
           IF (ldc) THEN
              CALL bsendc((it==0))
           ELSE
              CALL bsensi((it==0))
           END IF

!!!$   evtl   Rauhigkeitsmatrix belegen
           IF (l_bsmat) CALL bsmatm(it,l_bsmat)

        ELSE
!!!$   Felder zuruecksetzen
           sigma = sigma2

           sigmaa = sgmaa2

           IF (lrobust) wmatd = wmatd2

        END IF
        IF (itmax == 0) THEN ! only precalcs..
           errnr2 = 109
           PRINT*,'Only precalcs'
           EXIT
        END IF

!!!$   REGULARISIERUNG / STEP-LENGTH einstellen
        IF (.NOT.lstep) THEN
           IF (llam) THEN

!!!$   "Regularisierungsschleife" initialisieren und step-length zuruecksetzen
              llam = .FALSE.
              step = 1d0
           ELSE

!!!$   Regularisierungsindex hochzaehlen
              itr = itr+1
              IF (((((nrmsd.LT.rmsreg.AND.itr.LE.nlam).OR. &
                   (dlam.GT.1d0.AND.itr.LE.nlam)).AND.&
                   (.NOT.ldlamf.OR.dlalt.LE.1d0).AND.&
                   dabs(1d0-rmsreg/nrmsdm).GT.mqrms).OR.&
                   (rmsreg.EQ.0d0)).AND.&
                   (bdpar >= bdmin)) THEN
!!!$   Regularisierungsparameter bestimmen
                 IF (lsetup.OR.lsetip) THEN

!!!$   Kontrollausgabe
                    IF (llamf) THEN
                       lam = lamfix
                    ELSE
                       WRITE(*,'(a,i3,a,i3,a,t100,a)',ADVANCE='no')&
                            ACHAR(13)//' Iteration ',it,', ',itr,&
                            ' : Calculating 1st regularization parameter',''
                       WRITE(fprun,'(a,i3,a,i3,a)',ADVANCE='no')&
                            ' Iteration ',it,', ',itr,&
                            ' : Calculating 1st regularization parameter'
                       CALL blam0()
                       WRITE (*,'(a,G10.2)',ADVANCE='no')'lam_0:: ',lammax
                       WRITE (fprun,'(a,G10.2)')'lam_0 ',lammax
                       lam = lammax
!!!$   ak Model EGS2003, ERT2003                        call blam0()
!!!$   ak Model EGS2003, ERT2003                        lam = lammax
!!!$   ak                        lam = 1d4
                    END IF
                 ELSE
                    IF (llamf) THEN
                       lam = lamfix
                    ELSE
                       dlalt = dlam
                       IF (ldlami) THEN
                          ldlami = .FALSE.
                          alam   = dmax1(dabs(dlog(nrmsd/nrmsdm)),&
                               dlog(1d0+mqrms))
                          dlam   = fstart
                       ELSE
                          alam = dmax1(alam,dabs(dlog(nrmsd/nrmsdm)))
                          dlam = dlog(fstop)*&
                               SIGN(1d0,dlog(nrmsd/nrmsdm))+&
                               dlog(fstart/fstop)*&
                               dlog(nrmsd/nrmsdm)/alam
                          dlam = dexp(dlam)
                       END IF
                       lam = lam*dlam
                       IF (dlalt.GT.1d0.AND.dlam.LT.1d0) ldlamf=.TRUE.
!!!$   tst                        if (dlam.gt.1d0) lfstep=.true.
!!!$   ak Model EGS2003
                       IF (dlam.GT.1d0) lrobust=.FALSE.
                    END IF
                 END IF
              ELSE

!!!$   Regularisierungsparameter zuruecksetzen und step-length verkleinern
                 llam = .TRUE.
                 IF (llamf) THEN
                    lam = lamfix
                 ELSE
                    lam  = lam/dlam
                 END IF
                 IF (lfstep) THEN
                    lfstep = .FALSE.
                 ELSE
                    lstep = .TRUE.
                    step  = 5d-1
                 END IF
              END IF

!!!$   Ggf. Daten-RMS speichern
              IF (lsetup.OR.lsetip) THEN
                 lsetup = .FALSE.
                 lsetip = .FALSE.
              ELSE
                 IF (.NOT.lstep) rmsreg=nrmsd
              END IF
           END IF
        ELSE
           lstep = .FALSE.

!!!$   Parabolische Interpolation zur Bestimmung der optimalen step-length
           CALL parfit(rmsalt,nrmsd,rmsreg,nrmsdm,stpmin)

           IF (step.EQ.stpmin.AND.stpalt.EQ.stpmin)THEN

!!!$   Nach naechstem Modelling abbrechen
              stpalt = 0d0
           ELSE

!!!$   Step-length speichern
              stpalt = step
           END IF
           itr = itr + 1
        END IF
!!!$   Kontrollausgaben
        WRITE(*,'(a,i3,a,i3,a,t56,a)',ADVANCE='no')&
             ACHAR(13)//' Iteration ',it,', ',itr,&
             ' : Updating',''

        WRITE(fprun,*)' Iteration ',it,', ',itr,&
             ' : Updating'

!!!$ Modell parameter mit aktuellen Leitfaehigkeiten belegen
        CALL bpar
        IF (errnr.NE.0) GOTO 999

!!!$   UPDATE anbringen
        CALL update
        IF (errnr.NE.0) GOTO 999

!!!$ Leitfaehigkeiten mit verbessertem Modell belegen
        CALL bsigma
        IF (errnr.NE.0) GOTO 999

        IF (lverb) THEN
           CALL wout_up(kanal,it,itr)
           IF (errnr.NE.0) GOTO 999
        END IF
!!!$   Roughness bestimmen
        CALL brough

!!!$   Ggf. Referenzleitfaehigkeit bestimmen
        IF (lsr) CALL refsig()
!!!$   Neues Modelling
     END DO ! DO WHILE (.not. converged)

!!!$.................................................

!!!$ RESET FPI status variable to proceed with full COMPLEX calculus
     lip = .FALSE.
!!!$

!!!$   OUTPUT
     WRITE (*,'(a,t25,I4,t35,a,t100,a)')ACHAR(13)//&
          'MODEL ESTIMATE AFTER',it,'ITERATIONS',''
     CALL wout(kanal,dsigma,dvolt)
     IF (errnr.NE.0) GOTO 999

!!!$   Kontrollausgaben

     IF (errnr2.EQ.92) THEN

        WRITE(*,'(a22,a31)') ' Iteration terminated:',&
             ' Min. step-length for 2nd time.'

        WRITE(fprun,'(a22,a31)',err=999) ' Iteration terminated:',&
             ' Min. step-length for 2nd time.'
     ELSE IF (errnr2.EQ.80) THEN
        WRITE(*,'(a22,a10)') ' Iteration terminated:', &
             ' Min. RMS.'

        WRITE(fprun,'(a22,a10)',err=999) ' Iteration terminated:',&
             ' Min. RMS.'
     ELSE IF (errnr2.EQ.81) THEN
        WRITE(*,'(a22,a24)') ' Iteration terminated:',&
             ' Min. rel. RMS decrease.'

        WRITE(fprun,'(a22,a24)',err=999) ' Iteration terminated:',&
             ' Min. rel. RMS decrease.'
     ELSE IF (errnr2.EQ.79) THEN
        WRITE(*,'(a22,a19)') ' Iteration terminated:',&
             ' Max. # iterations.'

        WRITE(fprun,'(a22,a19)',err=999) ' Iteration terminated:',&
             ' Max. # iterations.'

     ELSE IF (errnr2.EQ.109) THEN
        WRITE(*,'(a)') ' Iteration terminated:'//&
             ' Min. model changes reached'

        WRITE(fprun,'(a)',err=999) ' Iteration terminated:'//&
             ' Min. model changes reached'

     END IF

!!!$   Run-time abfragen und ausgeben
     fetxt = ' CPU time: '
     CALL toc(c1,fetxt)
!!$     PRINT*,''
!!$     PRINT*,''
     WRITE(fprun,'(a)',err=999)TRIM(fetxt)

!!!$   Kontrolldateien schliessen
     CLOSE(fpinv)
     CLOSE(fpcjg)
     CLOSE(fpeps)

!!!$   Ggf. Summe der Sensitivitaeten aller Messungen ausgeben
     IF (lsens) THEN
        CALL BBSENS(kanal,dsens)
        IF (errnr.NE.0) GOTO 999
     END IF

     IF (lvario) THEN
        IF (lsens) CALL bvariogram_s ! calculate experimental variogram
        CALL bvariogram
     END IF

     IF (lcov1) CALL buncert (kanal,lamalt)

     CALL des_cjgmod(1,fetxt,errnr) ! call cjgmod destructor
     IF (errnr /= 0) GOTO 999

!!!$   'sens' und 'pot' freigeben
     IF (ldc) THEN
        fetxt = 'allocation sensdc'
        DEALLOCATE (sensdc)
        fetxt = 'allocation koptdc'
        DEALLOCATE (kpotdc)
     ELSE
        DEALLOCATE(sens,kpot)
     END IF

     fetxt = 'allocation smatm'
     IF (ALLOCATED (smatm)) DEALLOCATE (smatm)

     fetxt = 'allocation pot,pota,fak'
     IF (ALLOCATED (pot)) DEALLOCATE (pot,pota,fak)

     fetxt = 'allocation snr,sx,sy'
     IF (ALLOCATED (snr)) DEALLOCATE (snr,sx,sy)

     fetxt = 'allocation typ'
     IF (ALLOCATED (typ)) DEALLOCATE (typ,nelanz,selanz)

     fetxt = 'allocation nrel'
     IF (ALLOCATED (nrel)) DEALLOCATE (nrel,rnr)

     fetxt = 'allocation esp'
     IF (ALLOCATED (espx)) DEALLOCATE (espx,espy)

     fetxt = 'allocation kwn'
     IF (ALLOCATED (kwn)) DEALLOCATE (kwn)

     fetxt = 'allocation kwni'
     IF (ALLOCATED (kwnwi)) DEALLOCATE (kwnwi)

     fetxt = 'allocation elbg'
     IF (ALLOCATED (elbg)) DEALLOCATE (elbg,relbg,kg)

     fetxt = 'allocation enr'
     IF (ALLOCATED (enr)) DEALLOCATE (enr)

     fetxt = 'allocation mnr'
     IF (ALLOCATED (mnr)) DEALLOCATE (mnr)

     fetxt = 'allocation strnr,strom,volt,etc'
     IF (ALLOCATED (strnr)) DEALLOCATE (strnr,strom,volt,sigmaa,&
          kfak,wmatdr,wmatdp,vnr,dat,wmatd,wmatd2,sgmaa2,wdfak)

     fetxt = 'allocation par,dpar,dpar2'

     IF (ALLOCATED (par)) DEALLOCATE (par,dpar,dpar2)

     fetxt = 'allocation sigma'

     IF (ALLOCATED (sigma)) DEALLOCATE (sigma,sigma2)

     IF (ALLOCATED (d0)) DEALLOCATE (d0,fm0)
     IF (ALLOCATED (m0)) DEALLOCATE (m0)

     IF (ALLOCATED (rwddc)) DEALLOCATE (rwddc) 
     IF (ALLOCATED (rwndc)) DEALLOCATE (rwndc) 
     IF (ALLOCATED (rwd)) DEALLOCATE (rwd) 
     IF (ALLOCATED (rwn)) DEALLOCATE (rwn) 
     IF (ALLOCATED (rwdnr)) DEALLOCATE (rwdnr) 
     CLOSE(fprun)

!!!$   Ggf. weiteren Datensatz invertieren
  END DO

!!!$   'crtomo.cfg' schliessen
  CLOSE (fpcfg)

!!!$ crtomo.pid löschen
  fetxt = 'crtomo.pid'
  OPEN (fprun,FILE=TRIM(fetxt),STATUS='old',err=999)
  CLOSE (fprun,STATUS='delete')

!!!$ finished-string an inv.ctr anhängen
  fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
  OPEN(fpinv,file=TRIM(fetxt),status='old',POSITION='append',err=999)
  WRITE (fpinv,'(a)')'***finished***'
  CLOSE(fpinv)
  STOP '0'

!!!$.....................................................................

!!!$   Fehlermeldung
999 OPEN(fperr,file='error.dat',status='replace')
  errflag = 2
  CALL get_error(ftext,errnr,errflag,fetxt)
  WRITE(fperr,*) 'CRTomo PID ',pid,' exited abnormally'
  WRITE(fperr,'(a80,i3,i1)') fetxt,errnr,errflag
  WRITE(fperr,*)ftext
  CLOSE(fperr)

  STOP '-1'

END PROGRAM inv

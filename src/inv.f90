program inv

  ! Complex resistivity 2.5D inversion main program

  ! Used file id's  
  !  9 - error.dat   -> fperr
  ! 10 - run.ctr    -> fprun
  ! 11 - in-/output -> kanal
  ! 12 - crtomo.cfg -> fpcfg
  ! 13 - inv.ctr    -> fpinv
  ! 14 - cjg.ctr    -> fpcjg
  ! 15 - eps.ctr    -> fpeps
  !
  ! Andreas Kemna                                        02-May-1995
  ! Last change                                     	March 2014

  use alloci
  use tic_toc
  use femmod
  use datmod
  use invmod
  use cjgmod
  use sigmamod
  use electrmod
  use modelmod
  use elemmod
  use wavenmod
  use randbmod
  use konvmod
  use errmod
  use pathmod
  use bsmatm_mod
  use bmcm_mod
  use brough_mod
  use invhpmod
  use omp_lib
  use ompmod
  use get_ver
  use iso_fortran_env  

  implicit none
  character(256)        :: ftext
  integer               :: c1,i,count,mythreads,maxthreads
  real(prec)            :: lamalt
  logical               :: converged,l_bsmat
  integer               :: getpid,pid,myerr,info
  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  ! SETUP UND INPUT
  errnr = 1
  ! set file id's 
  fperr = 9 
  fprun = 10
  kanal = 11 
  fpcfg = 12 
  fpinv = 13 
  fpcjg = 14 
  fpeps = 15 

  ! Say hi
  call welcome()

  fetxt = 'crtomo.cfg'
  open(fpcfg,file=trim(fetxt),status='old',err=999)
  ! Switch: perform another inversion?
  lagain=.true. ! is set afterwards by user input file to false

  ! Set common parameters
  do while ( lagain ) ! loop over individual inversions, see 'lagain'

     errnr2 = 0
     ! Read all variables
     call rall(kanal,delem,delectr,dstrom,drandb,&
          dsigma,dvolt,dsens,dstart,dd0,dm0,dfm0,lagain)
     if (errnr.ne.0) goto 999

     !     call get_threads(nthreads,kwnanz)
     !    nthreads = 2
     !      CALL OMP_SET_NUM_THREADS ( NTHREADS )
     !   Element- und Randelementbeitraege sowie ggf. Konfigurationsfaktoren
     !   zur Berechnung der gemischten Randbedingung bestimmen
     call precal()

     if (errnr.ne.0) goto 999

     if (.not.lbeta) then
        lsr = .false.

        !   Ggf. Fehlermeldungen
        if (.not.lsink) then
           fetxt = 'no mixed boundary specify sink node'
           errnr = 102
        end if
        ! RM this is handeled in rall..
        if (swrtr.eq.0.and..not.lrandb2) then
           fetxt = ' '
           errnr = 103
        end if

        if (errnr.ne.0) goto 999
     else

        !   Ggf. Fehlermeldungen
        if (swrtr.eq.0) then
           fetxt = ' '
           errnr = 106
        end if  
     end if


     !   getting dynamic memory 
     errnr = 94
     ! physical model
     fetxt = 'allocation problem sigma'
     allocate (sigma(elanz),STAT=myerr)
     if (myerr /= 0) goto 999
     fetxt = 'allocation problem sigma2'
     allocate (sigma2(elanz),STAT=myerr)
     if (myerr /= 0) goto 999
     !  model parameters
     fetxt = 'allocation problem par'
     allocate (par(manz),STAT=myerr)
     if (myerr /= 0) goto 999
     fetxt = 'allocation problem dpar'
     allocate (dpar(manz),STAT=myerr)
     if (myerr /= 0) goto 999
     fetxt = 'allocation problem dpar2'
     allocate (dpar2(manz),STAT=myerr)
     if (myerr /= 0) goto 999
     fetxt = 'allocation problem pot'
     allocate(pot(sanz),STAT=myerr)
     if (myerr /= 0) goto 999
     fetxt = 'allocation problem pota'
     allocate (pota(sanz),STAT=myerr)
     if (myerr /= 0) goto 999
     ! now the big array are coming.. 
     fetxt = 'allocation problem fak'
     allocate (fak(sanz),STAT=myerr) ! fak for modeling
     if (myerr /= 0) goto 999
     if (ldc) then
        fetxt = 'allocation problem sensdc'
        allocate (sensdc(nanz,manz),STAT=myerr)
        if (myerr /= 0) goto 999
        fetxt = 'allocation problem kpotdc'
        allocate (kpotdc(sanz,eanz,kwnanz),STAT=myerr)
     else
        fetxt = 'allocation problem sens'
        allocate (sens(nanz,manz),STAT=myerr)
        if (myerr /= 0) goto 999
        fetxt = 'allocation problem kpot'
        allocate (kpot(sanz,eanz,kwnanz),STAT=myerr)
     end if
     if (myerr /= 0) goto 999

     ! get CG data storage of residuums and bvec, which is global
     call con_cjgmod (1,fetxt,myerr)
     if (myerr /= 0) goto 999

     ! >> RM ref model regu
     ! assign memory to global variables
     fetxt = 'allocation problem reference model'
     allocate (w_ref_re(manz),w_ref_im(manz),m_ref(manz),STAT=myerr)
     w_ref_re = 0d0;m_ref = cmplx(0d0);w_ref_im = 0d0
     if (myerr /= 0) goto 999
     ! << RM ref model regu


     !! INITIALIZE
     !   Startparameter setzen
     it     = 0;itr    = 0
     rmsalt = 0d0; lamalt = 1d0; bdpar = 1d0
     !     IF (lamnull_cri > 0d0) llamalt = lamnull_cri
     if (btest(llamf,0)) lamalt = lamfix
     betrms = 0d0; pharms = 0d0
     lsetup = .true.; lsetip = .false.; lfpi    = .false.
     llam   = .false.; ldlami = .true.; lstep  = .false.
     lfstep = .false.; l_bsmat = .true.
     step   = 1d0; stpalt = 1d0; alam   = 0d0

     !   Kontrolldateien oeffnen
     errnr = 1
     ! OPEN CONTRL FILES
     fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
     open(fpinv,file=trim(fetxt),status='replace',err=999)
     close(fpinv)
     fetxt = ramd(1:lnramd)//slash(1:1)//'run.ctr'
     open(fprun,file=trim(fetxt),status='replace',err=999)
     !  close(fprun) muss geoeffnet bleiben da sie staendig beschrieben wird
     fetxt = ramd(1:lnramd)//slash(1:1)//'cjg.ctr'
     open(fpcjg,file=trim(fetxt),status='replace',err=999)
     close(fpcjg)
     fetxt = ramd(1:lnramd)//slash(1:1)//'eps.ctr'
     open(fpeps,file=trim(fetxt),status='replace',err=999)

     !  SET ERRORS for all measurements and write control to fpeps
     ! >> RM
     if (ldc) then
        write (fpeps,'(a)')'1/eps_r      datum'
        write (fpeps,'(G12.3,2x,G14.5)')(sqrt(wmatdr(i)),&
             real(dat(i)),i=1,nanz)
     else

        write (fpeps,'(t5,a,t14,a,t27,a,t38,a,t50,a,t62,a,t71,a,t87,a)')'1/eps_r','1/eps_p',&
             '1/eps','eps_r','eps_p','eps','-log(|R|)', '-Phase (rad)'
        write (fpeps,'(3F10.1,2x,3G12.3,2G15.7)')&
             (sqrt(wmatdr(i)),sqrt(wmatdp(i)),sqrt(wmatd(i)),1/sqrt(wmatdr(i)),&
             1/sqrt(wmatdp(i)),1/sqrt(wmatd(i)),real(dat(i)),aimag(dat(i)),i=1,nanz)
     end if
     close(fpeps)
     errnr = 4


     ! set starting model 
     call bsigm0(kanal,dstart)
     if (errnr.ne.0) goto 999


     !   Kontrolldateien initialisieren
     !   diff-        call kont1(delem,delectr,dstrom,drandb)
     !   diff+<
     call kont1(delem,delectr,dstrom,drandb,dd0,dm0,dfm0,lagain)
     !   diff+>
     if (errnr.ne.0) goto 999

     !     write(6,"(a, i3)") " OpenMP max threads: ", OMP_GET_MAX_THREADS()
     !$OMP PARALLEL
     !     write(6,"(2(a,i3))") " OpenMP: N_threads = ",&
     !          OMP_GET_NUM_THREADS()," thread = ", OMP_GET_THREAD_NUM()
     !$OMP END PARALLEL

     !-------------
     write(*,*) '------------------------------------------'
     !   get current time
     call tic(c1)
     !.................................................
     converged = .false.

     do while (.not.converged) ! optimization loop

        !   Control output
        write(fprun,'(a,i3,a,i3,a)')' Iteration ',it,', ',itr,&
             ' : Calculating Potentials'
        !   MODELLING
        count = 0
        if (ldc) then
           fetxt = 'allocation problem adc'
           allocate (adc((mb+1)*sanz),STAT=myerr)
           if (myerr /= 0) goto 999
           fetxt = 'allocation problem hpotdc'
           allocate (hpotdc(sanz,eanz),STAT=myerr)
           if (myerr /= 0) goto 999
           allocate (bdc(sanz),STAT=myerr)
           fetxt = 'allocation problem adc'
           if (myerr /= 0) goto 999

           do k=1,kwnanz
              count = count + 1
              fetxt = 'DC-Calculation wavenumber'
              do l=1,eanz
                 if (lsr.or.lbeta.or.l.eq.1) then
                    !   Evtl calculation of analytical potentials
                    if (lsr) call potana(l,k,pota)
                    !   COMPilation of the linear system
                    fetxt = 'kompadc'
                    call kompadc(l,k,adc,bdc)
                    !   Evtl take Dirichlet boundary values into account
                    if (lrandb) call randdc(adc,bdc)
                    if (lrandb2) call randbdc2(adc,bdc)
                    !   Scale the linear system (preconditioning stores fak)
                    fetxt = 'scaldc'
                    call scaldc(adc,bdc,fak)
                    !   Cholesky-Factorization of the Matrix
                    fetxt = 'choldc'
                    call choldc(adc)
                 else
                    fetxt = 'kompbdc'
                    !   Modification of the current vector (Right Hand Side)
                    call kompbdc(l,bdc,fak)
                 end if
                 !   Solve linear system
                 fetxt = 'vredc'
                 call vredc(adc,bdc,pot)
                 !   Scale back the potentials, save them and 
                 !   eventually add the analytical response
                 do j=1,sanz
                    kpotdc(j,l,k) = real(pot(j)) * fak(j)
                    if (lsr) kpotdc(j,l,k) = kpotdc(j,l,k) + &
                         real(pota(j))
                    if (swrtr.eq.0) hpotdc(j,l) = kpotdc(j,l,k)
                 end do
              end do
           end do
        else
           fetxt = 'allocation problem a'
           allocate(a_mat_band(3*mb+1,sanz))
           allocate(a_mat_band_elec(3*mb+1,sanz))
           allocate (ipiv(sanz),STAT=myerr)
           if (myerr /= 0) goto 999
           fetxt = 'allocation problem hpot'
           allocate (hpot(sanz,eanz),STAT=myerr) 
           if (myerr /= 0) goto 999
           fetxt = 'allocation problem b'
           allocate (b(sanz),STAT=myerr)
           if (myerr /= 0) goto 999
           !   COMPLEX CASE
!!!$   POTENTIALWERTE BERECHNEN
           do k=1,kwnanz
!!!$   Kontrollausgabe
              count = count + 1

              if (swrtr.eq.0) then
                 write(*,'(a)')' Calculating Potentials'
              else
                 write (*,'(a,I3,a,I3,a,I3)',ADVANCE='no')&
                      achar(13)//' iteration',it,' update',itr,' wavenumver',count
              end if
              call pre_comp_ab(k,a_mat_band)

              do l=1,eanz
                 b = cmplx(0D0)
                 a_mat_band_elec = a_mat_band
                 b(enr(l)) = cmplx(1D0)
                 if (lsink) b(nsink) = cmplx(-1D0)
                 call comp_ab(k,a_mat_band_elec,l)
!!!!$   Ggf. Randbedingung beruecksichtigen
                 ! General Band matrix
                 call zgbsv(sanz,mb,mb, 1, a_mat_band_elec, 3*mb+1, ipiv, b, sanz,info )
                 if (info.ne.0) print*,'ZGBSV info:',info
!!!$   Potentialwerte zurueckskalieren und umspeichern sowie ggf.
!!!$   analytische Loesung addieren
                 do j=1,sanz
                    kpot(j,l,k) = b(j)
                    if (lsr) kpot(j,l,k) = kpot(j,l,k) + pota(j)
                    if (swrtr.eq.0) hpot(j,l) = kpot(j,l,k)
                 end do
              end do
           end do
        end if
        errnr = 0
        !   Ggf. Ruecktransformation der Potentialwerte
        if (swrtr.eq.1) then
           call rtrafo(errnr)
           if (errnr.ne.0) goto 999
        end if
        !   Spannungswerte berechnen
        call bvolti()
        if (errnr.ne.0) then
           !   reset model and data
           sigma = sigma2
           sigmaa = sgmaa2
           exit
        end if
        if (ldc) then
           deallocate(adc,hpotdc,bdc)
        else
           deallocate(hpot,b,ipiv,a_mat_band,a_mat_band_elec)
        end if
        if (lsetup.or.lsetip) then
           !   Ggf. background auf ratio-Daten "multiplizieren"
           if (lratio) then
              do j=1,nanz
                 dat(j) = dat(j) + sigmaa(j)
              end do
           end if
           !   Polaritaeten checken
           call chkpol(lsetup.or.lsetip)
        end if
        !   Daten-CHI berechnen
        call dmisft(lsetup.or.lsetip)
        if (errnr.ne.0) goto 999
        write(*,'(a,F8.2)') ' current chi**2:',nrmsd
        !   'nrmsd=0' ausschliessen
        if (nrmsd.lt.1d-12) nrmsd=nrmsdm*(1d0-mqrms)
        if (it == 0) then
           write(*,*) 'writing starting model'
           call wout(kanal,dsigma,dvolt)
        end if
        !   Kontrollvariablen ausgeben
        call kont2(lsetup.or.lsetip)
        if (errnr.ne.0) goto 999
        !   ABBRUCHBEDINGUNGEN
        if (llam.and..not.lstep) then
           !   Polaritaeten checken
           call chkpol(lsetup.or.lsetip)
           !   Wiederholt minimale step-length ?
           if (stpalt.eq.0d0) then
              errnr2 = 92
              fetxt = 'repeated step length'
           end if
           write (*,'(/a,G12.4,a/)')' convergence check (chi**2 (old/new)) ',&
                100.0*(1d0-rmsalt/nrmsd),'%'
           !   Keine Verbesserung des Daten-CHI ?
           if (abs(1d0-rmsalt/nrmsd).le.mqrms) then
              errnr2 = 81
              write (fetxt,*)'no CHI**2 improvement ',&
                   real(abs(1d0-rmsalt/nrmsd))
           end if
           !   Minimaler Daten-CHI erreicht ?
           !   tst            if (ABS(1d0-nrmsd/nrmsdm).le.mqrms) errnr2=80
           if (abs(1d0-nrmsd/nrmsdm).le.mqrms.and.ldlamf) then
              errnr2 = 80
              write (fetxt,*)'optimal chi**2 ',real(nrmsd),' reached'
           end if

           if (llam) then
              if (abs(1d0-nrmsd/rmsalt) < mqrms) errnr2 = 93
              if (nrmsd < 1d0 ) errnr2 = 94
              if (nrmsd > rmsalt) errnr2 = 95
           end if

           !   Maximale Anzahl an Iterationen ?
           if (it.ge.itmax) then
              errnr2 = 79
              write (fetxt,*)' reached max number of iterations ',itmax
           end if
           !   Minimal stepsize erreicht ?
           if (errnr2 == 0.and.bdpar <= bMIN) then
              errnr2 = 109
              write (fetxt,*)' Stepsize ',bdpar,' < Min stepsize ',bMIN
           end if

           !   Ggf. abbrechen oder "final phase improvement"
           if (errnr2.ne.0) then

              if (lfphai.and.errnr2.ne.79) then
                 print*,'CRI termination '//trim(fetxt),errnr2
                 !   ak
                 !   Widerstandsverteilung und modellierte Daten ausgeben
                 call wout(kanal,dsigma,dvolt)
                 if (errnr.ne.0) goto 999

                 !   Kontrollausgaben
                 write(*,'(/a/)')&
                      '**** Final phase improvement ****'

                 write(fprun,'(a24)',err=999)&
                      ' Final phase improvement'

                 fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
                 open(fpinv,file=trim(fetxt),status='old',&
                      POSITION='append',err=999)
                 write(fpinv,'(/a/)',err=999)&
                      '------------------------------------------------'//&
                      '------------------------------------------------'//&
                      '-----------------'
                 close(fpinv)

                 !   Wichtungsfeld umspeichern
                 wmatd = wmatdp
                 lam_cri = lamalt

                 write (*,'(/a,g12.4/)')'++ (FPI) setting phase error '//&
                      'and saving lam_cri: ',real(lam_cri)
                 write (fprun,'(/a,g12.4/)')'++ (FPI) setting phase error '//&
                      'and saving lam_cri: ',real(lam_cri)

                 lfpi    = .true.
                 lsetip = .true. ! 
                 lfphai = .false.
                 llam   = .false.
                 ldlami = .true.
                 lfstep = .true.
                 step   = 1d0
                 errnr2 = 0 ! reset convergence case "error"

                 !   ak
                 fetxt = 'cp -f inv.lastmod inv.lastmod_rho'
                 call SYSTEM (trim(fetxt))
                 if (lffhom) then
                    write(*,*)&
                         ' ******* Restarting phase model ********'
                    write(fprun,*)&
                         ' ******* Restarting phase model ********'
                    fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
                    open (fpinv,file=trim(fetxt),status='old',err=999,&
                         position='append')
                    write (fpinv,*)&
                         ' ******* Resetting phase model ********'
                    close (fpinv)
                    do j=1,elanz
                       sigma(j) = cmplx(&
                            cos(pha0/1d3)*abs(sigma(j)) ,&
                            -sin(pha0/1d3)*abs(sigma(j)) )
                    end do
                    lsetup = .true. ! ensure proper misfit and kont2 output
                    cycle       ! neues calc
                 else
                    lsetup = .false.
                 end if
                 !   ak

                 !   Daten-CHI berechnen
                 call dmisft(lsetip)
                 if (errnr.ne.0) goto 999
                 write (*,'(a,t45,a,t78,F14.4)') &
                      achar(13),'++ Phase data fit',nrmsd

                 !   Kontrollvariablen ausgeben
                 call kont2(lsetip)
                 if (errnr.ne.0) goto 999
              else
                 exit
              end if
           else
              ! >> RM
              lamalt = lam ! save the lambda of the previous iteration
              ! if, and only if the iterate was successful...
              ! << RM
              !   ak
              !   Widerstandsverteilung und modellierte Daten ausgeben
              write(*,*) 'writing model iterate',it                   
              call wout(kanal,dsigma,dvolt)
              if (errnr.ne.0) goto 999
           end if
        end if

        if ((llam.and..not.lstep).or.lsetup.or.lsetip) then
           !   Iterationsindex hochzaehlen
           it = it+1

           !   Parameter zuruecksetzen
           itr    = 0
           rmsreg = 0d0
           dlam   = 1d0
           dlalt  = 1d0
           ldlamf = .false.

           !   Daten-CHI speichern
           rmsalt = nrmsd

           !   Lambda speichen
           if (it>1) lamalt = lam ! mindestens einmal konvergiert

           !   Felder speichern
           sigma2 = sigma

           sgmaa2 = sigmaa

           if (lrobust) wmatd2 = wmatd

           !   Kontrollausgaben
           !           write (*,'(a)',ADVANCE='no')achar(13)//clear_line
           write (*,'(a,i3,a)') &
                achar(13)//' iteration',it,' sensitivities'

           write(fprun,'(a,i3,a,i3,a)')' Iteration ',it,', ',itr,&
                ' : Calculating Sensitivities'

           !   SENSITIVITAETEN berechnen
           if (ldc) then
              call bsendc((it==0))
           else
              call bsensi((it==0))
           end if

           !   evtl   Rauhigkeitsmatrix oder CmS-Matrix belegen
           if (l_bsmat) call bsmatm(it,l_bsmat)

        else
           if (lverb) then
              call wout_up(kanal,it,itr,.false.)
              if (errnr.ne.0) goto 999
           end if

           !   Felder zuruecksetzen
           sigma = sigma2

           sigmaa = sgmaa2

           if (lrobust) wmatd = wmatd2

        end if
        if (itmax == 0) then ! only precalcs..
           errnr2 = 109
           print*,'Only precalcs'
           exit
        end if

        if (btest(llamf,0)) then ! for fixed lambda we do not want any parabola fitting?
           lam = lamfix
           if (btest(llamf,1).and.it>1) lam = lamfix / (2d0 * real(it - 1))
           llam = .false. ! in order to calculate any update this needs to be false
           if (lsetup.or.lsetip) then
              lsetup = .false.
              lsetip = .false.
           end if

           !   REGULARISIERUNG / STEP-LENGTH einstellen
        else 
           if (.not.lstep) then

              if (llam) then
                 !   "Regularisierungsschleife" initialisieren und step-length zuruecksetzen
                 llam = .false.
                 step = 1d0
              else

                 !   Regularisierungsindex hochzaehlen
                 itr = itr+1
                 if ((((nrmsd.lt.rmsreg.and.itr.le.nlam).or. &
                      (dlam.gt.1d0.and.itr.le.nlam)).and.&
                      (.not.ldlamf.or.dlalt.le.1d0).and.&
                      (bdpar > bMIN).and.&
                      (abs(1d0-rmsreg/nrmsdm).gt.mqrms)).or.&
                      (rmsreg.eq.0d0)) then

                    if (rmsreg > 0d0) then
                       write (fprun,'(/a,G12.4,a)')'Chi increase:',&
                            100.0*(1d0-rmsalt/nrmsd),' %'
                       write (fprun,'(a,G12.4,a)')'Stepsize :',bdpar
                       write (fprun,'(a,G12.4/)')'nrmsd/rmsreg :',nrmsd/rmsreg
                    end if
                    !   Regularisierungsparameter bestimmen
                    if (lsetup.or.lsetip) then ! general initialization, lam0

                       !   Kontrollausgabe
                       write(*,'(a,i3,a,i3,a,t100,a)',ADVANCE='no')&
                            achar(13)//' iteration',it,' update',itr,&
                            ' lambda_0: '
                       write(fprun,'(a,i3,a,i3,a)',ADVANCE='no')&
                            ' Iteration ',it,' update',itr,&
                            ' : Calculating 1st regularization parameter'
                       call blam0()
                       !                       write (*,'(a,G10.2)',ADVANCE='no')'lam_0:: ',lammax
                       write (fprun,'(a,G10.2)')'lam_0 ',lammax
                       lam = lammax
                       lamalt = lammax
                       !   ak Model EGS2003, ERT2003                        call blam0()
                       !   ak Model EGS2003, ERT2003                        lam = lammax
                       !   ak                        lam = 1d4
                       !
                       ! GENERAL REMARK ON OUR REGULARISATION PARAMETER
                       ! Standard normal equations are
                       ! (A^h * Cd^-1 A + Cm^-1) dm = A_q^h * C_d^-1 * (d - f(m_q)) + C_m^-1 * ({m_q,(m_q-m_0)})
                       ! Where we identify lam as inverse a-priori variance:
                       ! Cm^-1 = \lam R^T * R.
                       ! If we solve the normal equation 
                       !  Cm * (A^h * Cd^-1 A + I) dm = Cm * (A_q^h * C_d^-1 * (d - f(m_q)) - ({m_q,(m_q-m_0)}))
                       ! instead, we use \lam as prior variance.
                       ! This also results in the fundamental problem on how to adjust the search direction which depends on the 
                       ! CHI decrease and thus should also be reciprocal.
                       !
                       ! nrmsdm and mqrms, fstart and fstop are set in rall.f90.
                       ! Defaults are
                       ! nrmsdm := 1d0; mqrms := 2d-2; fstart := 0.5; fstop := 0.9
                       !
                    else ! for lambda search we go here
                       dlalt = dlam
                       if (ldlami) then ! initialize search
                          ldlami = .false. 
                          alam   = max(abs(log(nrmsd/nrmsdm)),&
                               log(1.+mqrms))
                          ! alam = MAX(log(actual chi),log(1+0.02))
                          dlam   = fstart ! sets first dlam (0.5 default, which should be fine)
                       else
                          alam = max(alam,abs(log(nrmsd/nrmsdm))) 
                          ! CHI dependend partial fraction of lam variation
                          ! alam = MAX(alam,log(actual chi))
                          dlam = log(fstop)*&
                               sign(1._prec,log(nrmsd/nrmsdm))+&
                               log(fstart/fstop)*&
                               log(nrmsd/nrmsdm)/alam 
                          !
                          ! dlam = ln(0.9) * sign(1,ln(act chi)) + (ln(0.5/0.9)) * ln(act chi)/alam
                          ! this makes mostly the same and dlam = exp(ln(0.9) + ln(0.5/0.9)) = 0.5
                          ! This holds until the act chi value drops below 1. Than sign gives -1 and
                          ! also alam is now not identical to log(act chi) but log(1.02), which results
                          ! (for example act chi = 0,98)
                          ! dlam = exp(-ln(0.9) + ln(0.5/0.9) * ln(0.98)/log(1.02)) = 2
                          !
                          dlam = exp(dlam)
                       end if
                       lam = lam*dlam
                       if (dlalt.gt.1d0.and.dlam.lt.1d0) ldlamf=.true.
                       !   tst                        if (dlam.gt.1d0) lfstep=.true.
                       !   ak Model EGS2003
                       if (dlam.gt.1d0) lrobust=.false.

                    end if
                 else

                    !   Regularisierungsparameter zuruecksetzen und step-length verkleinern
                    ! if no Chi decrease found
                    llam = .true.
                    lam  = lam/dlam

                    if (lfstep) then
                       lfstep = .false.
                    else
                       lstep = .true.
                       step  = 5d-1
                    end if
                 end if

                 !   Ggf. Daten-CHI speichern
                 if (lsetup.or.lsetip) then
                    lsetup = .false.
                    lsetip = .false.
                 else
                    if (.not.lstep) rmsreg=nrmsd
                 end if
              end if
           else
              lstep = .false.

              !   Parabolische Interpolation zur Bestimmung der optimalen step-length
              call parfit(rmsalt,nrmsd,rmsreg,nrmsdm,stpmin)

              if (step.eq.stpmin.and.stpalt.eq.stpmin)then

                 !   Nach naechstem Modelling abbrechen
                 stpalt = 0d0
              else

                 !   Step-length speichern
                 stpalt = step
              end if

           end if

        end if
        !   Kontrollausgaben
        !        write (*,'(a)',ADVANCE='no')achar(13)//clear_line
        !        write(*,'(a,i3,a,i3,a)',ADVANCE='no')&
        !             achar(13)//' Iteration ',it,', ',itr,&
        !             ' : Updating'

        write(fprun,*)' Iteration ',it,', ',itr,&
             ' : Updating'

        ! Modell parameter mit aktuellen Leitfaehigkeiten belegen

        call bpar
        if (errnr.ne.0) goto 999

        !   UPDATE anbringen
        call update
        if (errnr.ne.0) goto 999

        ! Leitfaehigkeiten mit verbessertem Modell belegen
        call bsigma
        if (errnr.ne.0) goto 999

        if (lverb) then
           call wout_up(kanal,it,itr,.true.)
           if (errnr.ne.0) goto 999
        end if
        !   Roughness bestimmen
        call brough

        !   Ggf. Referenzleitfaehigkeit bestimmen
        ! >>> RM
        !! $THIS has now a general meaning with the 
        ! mixed boundary, since we like to set sigma0
        ! as reference sigma as "mean" boundary value
        ! <<< RM
        if (lbeta) call refsig()

        if (btest(llamf,0)) then
           llam = .true.
        end if

        !   Neues Modelling
     end do ! DO WHILE (.not. converged)

     ! RESET FPI status variable to proceed with full COMPLEX calculus
     lfpi = .false.
     !
     !.................................................

     !   OUTPUT
     write (*,*)&
          'Model estimate after',it,'iterations'


     !   Kontrollausgaben
     ! errnr2 is more a status variable than a pure error number
     ! it shows which case was the reason for termination of the algorithm
     select case (errnr2)
     case (95)

        write(*,'(a22,a31)') ' Iteration terminated:',&
             ' no CHI decrease'
        write(fprun,'(a22,a31)',err=999) ' Iteration terminated:',&
             ' no CHI decrease'
        ! reset model to previous state
        if (btest(llamf,0)) then
           write (*,'(a)',ADVANCE='no')'Taking model state of previous iteration... '
           sigma = sigma2
           sigmaa = sgmaa2
           nrmsd = rmsalt
           print*,'CHI = ',real(nrmsd)
        end if

     case (94)

        write(*,'(a22,a31)') ' Iteration terminated:',&
             ' CHI < 1'
        write(fprun,'(a22,a31)',err=999) ' Iteration terminated:',&
             ' CHI < 1'

     case (93)

        write(*,'(a22,a25)') ' Iteration terminated:',&
             ' CHI decrease sufficient'
        write(fprun,'(a22,a25)',err=999) ' Iteration terminated:',&
             ' CHI decrease sufficient'

     case (92)

        write(*,'(a22,a31)') ' Iteration terminated:',&
             ' Min. step-length for 2nd time.'

        write(fprun,'(a22,a31)',err=999) ' Iteration terminated:',&
             ' Min. step-length for 2nd time.'
     case (80)
        write(*,'(a22,a10)') ' Iteration terminated:', &
             ' Min. CHI.'

        write(fprun,'(a22,a10)',err=999) ' Iteration terminated:',&
             ' Min. CHI.'
     case (81)
        write(*,'(a22,a24)') ' Iteration terminated:',&
             ' Min. rel. CHI decrease.'

        write(fprun,'(a22,a24)',err=999) ' Iteration terminated:',&
             ' Min. rel. CHI decrease.'
     case (79)
        write(*,'(a22,a19)') ' Iteration terminated:',&
             ' Max. # iterations.'

        write(fprun,'(a22,a19)',err=999) ' Iteration terminated:',&
             ' Max. # iterations.'

     case (109)
        write(*,'(a)') ' Iteration terminated:'//&
             ' Min. model changes reached'

        write(fprun,'(a)',err=999) ' Iteration terminated:'//&
             ' Min. model changes reached'

     end select

     call wout(kanal,dsigma,dvolt)
     if (errnr.ne.0 .and..not. errnr == 82) goto 999

     !   Kontrollausgaben

     !   Run-time abfragen und ausgeben
     fetxt = 'CPU time: '
     call toc(c1,fetxt)
     print*,''
     print*,''
     write(fprun,'(a)',err=999)trim(fetxt)

     !   Kontrolldateien schliessen
     close(fpinv)
     close(fpcjg)
     close(fpeps)


     !   Ggf. Summe der Sensitivitaeten aller Messungen ausgeben
     if (lsens) then
        call BBSENS(kanal,dsens)
        if (errnr.ne.0) goto 999
     end if

     if (lvario) then
        !        IF (lsens) CALL bvariogram_s ! calculate experimental variogram
        call bvariogram
     end if

     if (lcov1) call buncert (kanal,lamalt)

     call des_cjgmod(1,fetxt,myerr) ! call cjgmod destructor
     if (myerr /= 0) goto 999

     !   'sens' und 'pot' freigeben
     if (ldc) then
        fetxt = 'allocation sensdc'
        deallocate (sensdc)
        fetxt = 'allocation koptdc'
        deallocate (kpotdc)
     else
        deallocate(sens,kpot)
     end if

     fetxt = 'allocation smatm'
     if (allocated (smatm)) deallocate (smatm)

     fetxt = 'allocation pot,pota,fak'
     if (allocated (pot)) deallocate (pot,pota,fak)

     fetxt = 'allocation snr,sx,sy'
     if (allocated (snr)) deallocate (snr,sx,sy)

     fetxt = 'allocation typ'
     if (allocated (typ)) deallocate (typ,nelanz,selanz)

     fetxt = 'allocation nrel'
     if (allocated (nrel)) deallocate (nrel,rnr)

     fetxt = 'allocation esp'
     if (allocated (espx)) deallocate (espx,espy)

     fetxt = 'allocation kwn'
     if (allocated (kwn)) deallocate (kwn)

     fetxt = 'allocation kwni'
     if (allocated (kwnwi)) deallocate (kwnwi)

     fetxt = 'allocation elbg'
     if (allocated (elbg)) deallocate (elbg,relbg,kg)

     fetxt = 'allocation enr'
     if (allocated (enr)) deallocate (enr)

     fetxt = 'allocation mnr'
     if (allocated (mnr)) deallocate (mnr)

     fetxt = 'allocation strnr,strom,volt,etc'
     if (allocated (strnr)) deallocate (strnr,strom,volt,sigmaa,&
          kfak,wmatdr,wmatdp,vnr,dat,wmatd,wmatd2,sgmaa2,wdfak,&
          wmatd_cri) !!! these are allocated in rdati!!!

     fetxt = 'allocation par,dpar,dpar2'

     if (allocated (par)) deallocate (par,dpar,dpar2)

     fetxt = 'allocation sigma'

     if (allocated (sigma)) deallocate (sigma,sigma2)

     if (allocated (d0)) deallocate (d0,fm0)
     if (allocated (m0)) deallocate (m0)

     if (allocated (rwddc)) deallocate (rwddc) 
     if (allocated (rwndc)) deallocate (rwndc) 
     if (allocated (rwd)) deallocate (rwd) 
     if (allocated (rwn)) deallocate (rwn) 
     if (allocated (rwdnr)) deallocate (rwdnr) 

     if (allocated (w_ref_re)) deallocate (w_ref_re,w_ref_im,m_ref)

     close(fprun)

     !   Ggf. weiteren Datensatz invertieren
  end do

  !   'crtomo.cfg' schliessen
  close (fpcfg)

  ! crtomo.pid löschen
  fetxt = 'crtomo.pid'
  open (fprun,FILE=trim(fetxt),STATUS='old',err=999)
  close (fprun,STATUS='delete')

  ! finished-string an inv.ctr anhängen
  fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
  open(fpinv,file=trim(fetxt),status='old',POSITION='append',err=999)
  write (fpinv,'(a)')'***finished***'
  close(fpinv)
  stop '0'

  !.....................................................................

  !   Fehlermeldung
999 open(fperr,file='error.dat',status='replace')
  errflag = 2
  call get_error(ftext,errnr,errflag,fetxt)
  write(fperr,*) 'CRTomo PID ',pid,' exited abnormally'
  write(fperr,'(a80,i3,i1)') fetxt,errnr,errflag
  write(fperr,*)ftext
  close(fperr)

  stop '-1'

end program inv

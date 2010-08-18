      subroutine rall(kanal,delem,delectr,dstrom,drandb,
c     diff-     1                  dsigma,dvolt,dsens,dstart,lsens,lagain)
c     diff+<
     1     dsigma,dvolt,dsens,dstart,dd0,dm0,dfm0,lagain)
c     diff+>
      
c     Unterprogramm zum Einlesen der benoetigten Variablen.

c     Andreas Kemna                                            01-Mar-1995
c     Letzte Aenderung   20-Aug-2007
c     
c.....................................................................

      USE make_noise
      USE variomodel
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
      USE errmod
      USE konvmod
      USE pathmod

      IMPLICIT none


c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Kanalnummer
      integer         * 4     kanal

c     Dateinamen
      character       * 80    delem,
     1     delectr,
     1     dstrom,
     1     dsigma,
     1     dvolt,
     1     dsens,
     1     dstart,
c     diff+<
     1     dd0,
     1     dm0,
     1     dfm0,
c     diff+>
     1     drandb

c     Schalter ob weiterer Datensatz invertiert werden soll
      logical         * 4     lagain
      logical         * 4     lsto
c     check whether the file format is crtomo konform or not..
      logical           ::    crtf

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Indexvariable
      integer         * 4     i,iregus,ifp1

c     Pi
      real            * 8     pi

c     diff+<
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: dum,dum2
      REAL(KIND(0D0))                            :: dum3
      INTEGER(KIND = 4),DIMENSION(:),ALLOCATABLE :: idum,ic,ip
      integer         * 4     nanz0,j,j0
c     diff+>

c     ak Inga
      integer         * 4     elec1,elec2,
     1     elec3,elec4
      character(120) :: buff
      logical        :: exi
c.....................................................................

      pi = dacos(-1d0)

c     'crtomo.cfg' EINLESEN
      fetxt = 'crtomo.cfg'
      errnr = 3
c#################DEFAULTS ########################
c###  switches 
c     ro        lsr    = .false.
c     ro        lpol   = .true.
c     ro        lfphai = .true.
c     ro        lrho0  = .false.
c     akERT2003
c     ak        ldc    = .true.
c     ak        lsr    = .false.
c     ak        lpol   = .false.
c     Sonstiges
c     diff+<
      lsr    = .false.
      lpol   = .false.
      lindiv = .false.
c     ak
c     'dstart'
      lstart = .false.
      dstart = ' '
c     ak        lstart = .true.
c     ak        dstart = '..\..\strasbrg\9610\plane45\mod\rho0.dat'
      ltri   = 0
      lsto = .false.            !default--
c     "Force negative phase" ?
c     sandra        lphi0 = .true.
      lphi0 = .FALSE.
c     ak        lphi0 = .false.
c     "ratio-dataset" ?
      lratio = .false.
c     ak        lratio = .true.
      llamf = .FALSE.
c     final phase improvement setzt phase zueruck auf homogenes modell
      lffhom = .FALSE.
c     Daten Rauschen vom Fehlermodell entkoppeln ?
      lnse2 = .FALSE.
c######values..
c     FIXED PARAMETER
c     Slash
      slash = '/'
c     Minimale "L1-ratio" (Grenze der "robust inversion")
      l1min = 1d0
c     ak        l1min = 1.2d0
      nrmsdm = 1d0
c     Art der Ruecktransformation
c     ak        swrtr = 1
c     Minimaler Quotient zweier aufeinanderfolgender Daten-RMS-Werte
c     ak Default
c     ak        mqrms = 1d-2
c     ak ERT 2002-2003 synth
      mqrms = 2d-2
c     ak        mqrms = 2d-2
c     ak Tank
c     ak        mqrms = 2d-2
c     ak MMAJ
c     ak        mqrms = 5d-2
c     CG-Epsilon
      eps = 1d-4
c     Mindest-step-length
      stpmin = 1d-3
c     Minimale stepsize (bdpar)
      bdmin = 0.0d0
c     Regularisierungsparameter
c     ak Default
      nlam   = 30
c     ak Default
      fstart = 0.5d0
      fstop  = 0.9d0
      lamfix = 0.0D0
c     ak MMAJ
c     ak        fstart = 0.2d0
c     ak        fstop  = 0.8d0
c     ak Strasbrg/Werne/Grimberg
c     ak        fstart = 0.5d0
c     ak        fstop  = 0.8d0
      iseedpri = 0; modl_stdn = 0.; iseed = 1;
      mswitch = 0
      iregus = 0
c#########################################################
c     Read in input values..
      fetxt = 'rall -> grid file'
      read(fpcfg,*,end=1001,err=98) mswitch
 98   read(fpcfg,'(a80)',end=1001,err=999) delem
      fetxt = 'rall -> electrode file'
      read(fpcfg,'(a80)',end=1001,err=999) delectr
      fetxt = 'rall -> spannungs file'
      read(fpcfg,'(a80)',end=1001,err=999) dstrom
      fetxt = 'rall -> Inversionsverzeichnis'
      read(fpcfg,'(a80)',end=1001,err=999) ramd
      INQUIRE (FILE=TRIM(ramd),EXIST= exi)
      IF (.NOT.exi) CALL SYSTEM ('mkdir '//TRIM(ramd))
c     diff+<
      fetxt = 'rall -> Differenz inversion'
      read(fpcfg,*,end=1001,err=999) ldiff
      fetxt = 'rall -> Diff. Messspannung'
      read(fpcfg,'(a80)',end=1001,err=999) dd0
      fetxt = 'rall -> Diff. Modell (auch prior)'
      read(fpcfg,'(a80)',end=1001,err=999) dm0
      fetxt = 'rall -> Diff. Modellierungen'
      read(fpcfg,'(a80)',end=1001,err=999) dfm0

      IF (dm0 /= '') THEN
         INQUIRE(FILE=TRIM(dm0),EXIST=lstart) ! prior model ?
         IF (lstart) THEN       ! set the starting model
            dstart = dm0
            PRINT*,'reading prior:',ACHAR(9)//TRIM(dm0)
         ELSE
            PRINT*,'omitting prior:',ACHAR(9)//TRIM(dm0)
            dm0 = ''
         END IF
      END IF
      IF (lstart.AND.ldiff.AND.((dd0 == ''.AND.dfm0 == ''))) THEN
         PRINT*,'Reference model regularization!'
         lprior = .TRUE.        ! reference model regu only if there is no
         ldiff = .FALSE.        ! time difference inversion
      END IF
c     diff+>
      fetxt = 'trying noise model seed'
      read(fpcfg,*,end=1001,err=99) iseedpri,modl_stdn
!     hier landet man nur, wenn man iseed und modl_stdn angenommen hat
      lnse2 = .NOT.lprior       ! kein prior?
!     Daten Rauschen unabhängig vom Fehlermodell?
      lnsepri = lprior          ! if we have seed and std we assume to add noise to prior
 99   fetxt = 'rall -> Gitter nx'
      read(fpcfg,*,end=1001,err=999) nx
      fetxt = 'rall -> Gitter nz'
      read(fpcfg,*,end=1001,err=999) nz
      fetxt = 'rall -> Anisotropie /x'
      read(fpcfg,*,end=1001,err=999) alfx
      fetxt = 'rall -> Anistotropie /y'
      read(fpcfg,*,end=1001,err=999) alfz
      fetxt = 'rall -> Maximale Iterationen'
      read(fpcfg,*,end=1001,err=999) itmax
c     ak        read(fpcfg,*,end=1001,err=999) nrmsdm
      fetxt = 'rall -> DC/IP Inversion'
      read(fpcfg,*,end=1001,err=999) ldc
c     ak        read(fpcfg,*,end=1001,err=999) lsr
      fetxt = 'rall -> Robuste Inversion'
      read(fpcfg,*,end=1001,err=999) lrobust
c     ak        read(fpcfg,*,end=1001,err=999) lpol
      fetxt = 'rall -> Finale Phasen Inversion'
      read(fpcfg,*,end=1001,err=999) lfphai
c     ak        read(fpcfg,*,end=1001,err=999) lindiv
      fetxt = 'rall -> Relativer Fehler Widerstand'
      read(fpcfg,*,end=1001,err=999) stabw0
      fetxt = 'rall -> Absoluter Fehler Widerstand'
      read(fpcfg,*,end=1001,err=999) stabm0
      fetxt = 'rall -> Phasenfehlerparameter A1'
      read(fpcfg,*,end=1001,err=999) stabpA1
      fetxt = 'rall -> Phasenfehlerparameter B'
      read(fpcfg,*,end=1001,err=999) stabpB
      fetxt = 'rall -> Relative Fehler Phasen'
      read(fpcfg,*,end=1001,err=999) stabpA2
      fetxt = 'rall -> Absoluter Fehler Phasen (mRad)'
      read(fpcfg,*,end=1001,err=999) stabp0
      fetxt = 'rall -> Homogenes Startmodell?'
      read(fpcfg,*,end=1001,err=999) lrho0
      fetxt = 'rall -> rho_0'
      read(fpcfg,*,end=1001,err=999) bet0
      fetxt = 'rall -> phase_0'
      read(fpcfg,*,end=1001,err=999) pha0
      fetxt = 'rall -> Noch eine Inversion'
      read(fpcfg,*,end=1001,err=999) lagain
      fetxt = 'rall -> 2D oder 2.5D ?'
      read(fpcfg,*,end=1001,err=999) swrtr
      fetxt = 'rall -> weitere Quelle?'
      read(fpcfg,*,end=1001,err=999) lsink
      fetxt = 'rall -> Nummer der Quelle'
      read(fpcfg,*,end=1001,err=999) nsink
      fetxt = 'rall -> Randbedingungen ?'
      read(fpcfg,*,end=1001,err=999) lrandb2
      fetxt = 'rall -> Datei mit Randwerten'
      read(fpcfg,'(a80)',end=1001,err=999) drandb
      fetxt = 'triangularization switch'
      read(fpcfg,'(I2)',end=100,err=100) ltri

      IF (ltri >= 20) THEN
         llamf = .TRUE.
         READ(fpcfg,*,end=104,err=104) lamfix
	 GOTO 105
 104     lamfix = 1.0           ! default value for MGS
         BACKSPACE(fpcfg)

 105     PRINT*,'Fixing Lambda =', lamfix
         ltri = ltri - 20
      END IF

      lsto = (ltri == 15)
      
      GOTO 101

 100  BACKSPACE (fpcfg)


 101  IF (lsto) PRINT*,'Stochastische Regularisierung'
      
      IF (ltri > 4 .AND. ltri < 15) THEN
         READ(fpcfg,*,end=102,err=102) betamgs
	 GOTO 103
 102     betamgs = 0.1          ! default value for MGS
         BACKSPACE (fpcfg)

 103     PRINT*,'Regularisation with support stabilizer beta =',
     1        betamgs
      END IF      

      IF (itmax == 0) PRINT*,
     1     ' ####### Only precalcs, itmax==0 ###########'

c     check if the final phase should start with homogenous model      
      lffhom = (stabp0 < 0)
      IF (lffhom) stabp0 = -stabp0
      
      lnse = ( stabw0 < 0 )     ! couple error and noise model
      IF ( lnse ) THEN
         stabw0 = -stabw0
         IF (lnse2) print*,'overriding seperate noise model'
         lnse2 = .FALSE.        ! overrides the lnse2 switch
c     copy error model into noise model
         nstabw0 = stabw0
         nstabm0 = stabm0
         nstabpA1 = stabpA1
         nstabpA2 = stabpA2
         nstabp0 = stabp0

         READ(fpcfg,*,end=106,err=106) iseed
	 GOTO 107
 106     iseed = 1              ! default value for PRS
         BACKSPACE(fpcfg)
         WRITE (*,'(a)')' Rauschen '//
     1        'Gekoppelt an Fehlermodell '
      ELSE
c     check if there is at least crt.noisemod containig noise info
         IF (iseedpri == 0) iseedpri = 1
         fetxt = 'crt.noisemod'
         INQUIRE(FILE=TRIM(fetxt),EXIST=lnse2)
      END IF
 107  IF (lnse2) THEN

         iseed = iseedpri
         WRITE (*,'(a,I7)',ADVANCE='no')
     1        'Entkoppeltes Daten Rauschen:: seed:',iseed
         
         nstabw0 = modl_stdn
         
         fetxt = 'get noise model from crt.noisemod'
         CALL get_noisemodel(nstabw0,nstabm0,nstabpA1,
     1        nstabpB,nstabpA2,nstabp0,errnr)
         IF (errnr /= 0) GOTO 999

         modl_stdn = 0.
         iseedpri = 0
         
         lnse = .TRUE.          ! add noise

      END IF

      IF (lnse) THEN 
         fetxt = 'write out noise model'
         CALL write_noisemodel(nstabw0,nstabm0,
     1        nstabpA1,nstabpB,nstabpA2,nstabp0,errnr)
         IF (errnr /= 0) GOTO 999
      ELSE
         PRINT*,'No Data noise!!'
      END IF

      IF ((nx<=0.OR.nz<=0).AND.ltri==0) ltri=1 ! at least L1-smoothness

c     Ggf. Fehlermeldungen
      if (ltri==0.AND.(nx.lt.2.or.nz.lt.2)) then
         fetxt = ' '
         errnr = 89
         goto 999
c$$$  else if (alfx.le.0d0.or.alfz.le.0d0) then
c$$$  fetxt = ' '
c$$$  errnr = 96
c$$$  goto 999
      else if (itmax<0.or.itmax.ge.100) then
         fetxt = ' '
         errnr = 61
         goto 999
      else if (nrmsdm.lt.1d-12) then
         fetxt = ' '
         errnr = 62
         goto 999
c     else if (nlam.lt.0) then
c     fetxt = ' '
c     errnr = 83
c     goto 999
c     else if (fstart.gt.1d0.or.fstop.gt.1d0.or.
c     1           fstart.le.0d0.or.fstop.le.0d0.or.
c     1           fstart.gt.fstop) then
c     fetxt = ' '
c     errnr = 98
c     goto 999
      else if (stabw0.le.0d0.or.stabm0.lt.0d0) then
         fetxt = ' '
         errnr = 104
         goto 999
      else if (.not.ldc.and.lfphai.and.
     1        ((stabp0.lt.0d0.or.stabpA2.lt.0d0).OR.
     1        ((stabp0 == 0d0).and.(stabpA2 == 0d0)))) then
         fetxt = ' '
         errnr = 105
         goto 999
      else if (lrho0.and.(bet0.le.0d0.or.
     1        (.not.ldc.and.dabs(pha0).gt.1d3*pi))) then
         fetxt = ' '
         errnr = 91
         goto 999
c     else if (mqrms.lt.0d0.or.mqrms.ge.1d0) then
c     fetxt = ' '
c     errnr = 64
c     goto 999
c     else if (lrobust.and.l1min.lt.1d0) then
c     fetxt = ' '
c     errnr = 90
c     goto 999
      end if

      lelerr = .NOT.lfphai.AND..NOT.ldc ! complex inversion only

c     (mswitch) Mega switch testing..
      lsens = BTEST(mswitch,0)  ! +1 ueberdeckung schreiben
      lcov1 = BTEST(mswitch,1)  ! +2 posterior modell covariance matrix 1
      lres  = BTEST(mswitch,2)  ! +4 rsolution matrix berechnen
      lcov2 = BTEST(mswitch,3)  ! +8 posterior modell covariance matrix 2
      
      lgauss = BTEST (mswitch,4) ! +16 solve ols with Gauss elemination
      
      lelerr = BTEST (mswitch,5).OR.lelerr ! +32 overwrites previous lelerr
      
      lsytop = .NOT.BTEST (mswitch,8) ! +256 disables sy top check of 
!     no flow boundary electrodes for enhanced beta calculation (bsytop). 
!     This is useful for including topographical effects and should be used
      
      lverb = BTEST (mswitch,10) ! +1024 Verbose output CG, daten, bnachbar..
      
      IF (lverb) WRITE(*,'(/a/)')' #  ## VERBOSE ## #'

      lres = (lres.or.lcov2)    ! compute mcm2 on top of resolution
      lcov1 = (lres.or.lcov1)   ! compute resolution by taking mcm1
c     
      lsens = .TRUE.            ! default immer coverages schreiben..
c     
      if (lratio) then
         lrho0  = .true.
         lstart = .false.
         lphi0  = .false.
         lpol   = .false.
      end if
c     diff-        if (lstart) lrho0=.false.

      if (lstart.or.ldiff) lrho0=.false.
c     ak
      if (ldiff) then
         ldc  = .true.
         lpol = .false.
      end if
c     diff+>
c     ak        if (ldc.or.stabp0.ge.stabw0) lfphai=.false.
      if (ldc) lfphai=.false.
c     Dateien
      lnramd = index(ramd,' ')-1
      dsigma = ramd(1:lnramd)//slash(1:1)//'rho.dat'
      dvolt  = ramd(1:lnramd)//slash(1:1)//'volt.dat'
      dsens  = ramd(1:lnramd)//slash(1:1)//'coverage.mag'

c     Elementeinteilung einlesen
      WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//'reading grid'
      call relem(kanal,delem)
      if (errnr.ne.0) goto 999

      IF (ltri/=0) THEN
         manz = elanz           ! wichtig an dieser stelle..
         CALL bnachbar          ! blegt nachbar
         CALL besp_elem
         lvario = .TRUE.
      ELSE
c     Modelleinteilung gemaess Elementeinteilung belegen
         manz = nx*nz           ! nur für strukturierte gitter
      END IF
      IF (lstart .OR. ldiff .OR. lprior) THEN
         ALLOCATE (m0(manz),stat=errnr)
         IF (errnr /= 0) THEN
            fetxt = 'Error memory allocation m0'
            errnr = 94
            goto 999
         END IF
      END IF
      lvario = lvario.OR.       ! if already set or
     1     (itmax == 0).AND.(lstart.OR.lprior) ! analyse any prior

      IF (lvario) CALL set_vario (nx,alfx,alfz,esp_mit,esp_med) ! nx is than
!     the variogram and covariance function type, see variomodel.f90

      if (manz.ne.elanz) then
         fetxt = 'manz /= elanz .. is not implemented yet'
         errnr = 50
         goto 999
      end if
!     !$ get memory for mnr..
      ALLOCATE (mnr(elanz),stat=errnr)
      IF (errnr /= 0) THEN
         fetxt = 'Error memory allocation mnr failed'
         errnr = 94
         goto 999
      END IF
!     !$ set mnr.. this may be altered if we have zonal approach..
      do i=1,elanz
         mnr(i) = i
      end do
!     !$ all the model mapping is than beased on a new numbering..
!     !$ zonal approach ?

c     Maximale Anzahl an CG-steps setzen
c     ak        ncgmax = manz
      ncgmax = manz / 2         ! useful for small scale model variations
!     for normal smooth and damping we usually need fewer CG iterations;
!     because the model variations are of bigger scale size
      IF (ltri < 5) ncgmax = ncgmax / 10
c     Elektrodenverteilung und Daten einlesen sowie Wellenzahlwerte
c     bestimmen

      call relectr(kanal,delectr)
      if (errnr.ne.0) goto 999

      IF (lsink) THEN
         WRITE(*,'(/A,I5,2F12.3/)')'Fictious sink @ node ',
     1        nsink,sx(snr(nsink)),sy(snr(nsink))
         WRITE(fpinv,'(A,I5,2F12.3)')'Fictious sink @ node ',
     1        nsink,sx(snr(nsink)),sy(snr(nsink))
      END IF

      call rdati (kanal,dstrom)

      if (errnr.ne.0) goto 999

      if (swrtr.eq.0) then
         lsr    = .false.
         kwnanz = 1
         kwn(1) = 0d0
         do i=1,typanz
            IF (typ(i) == 11) THEN
               PRINT*,'hier am besten aussteigen '//
     1              'da es im 2D keine gemischten RB gibt'
               fetxt = 'hier am besten aussteigen '//
     1              'da es im 2D keine gemischten RB gibt'
               errnr = 110
               GOTO 999
            END IF
         END DO
      else
         call rwaven()
         if (errnr.ne.0) goto 999
      end if

c     read boundary values
      if (lrandb2) then
         call rrandb(kanal,drandb)
         if (errnr.ne.0) goto 999
      end if
c     diff+<
      if (ldiff) then
         ALLOCATE (d0(nanz),fm0(nanz),stat=errnr)
         IF (errnr /= 0) THEN
            fetxt = 'Error memory allocation diff data '
            errnr = 94
            goto 999
         END IF
         open(kanal,file=dd0,status='old')
         read(kanal,*) nanz0
         read(kanal,*,err=999) elec1
         BACKSPACE(kanal)

         elec3=elec1-10000      ! are we still positive?
         crtf=(elec3 > 0)       ! crtomo konform?
         
         ALLOCATE (dum(nanz0),dum2(nanz0),idum(nanz0),
     1        ic(nanz0),ip(nanz0),stat=errnr)
         IF (errnr /= 0) THEN
            fetxt = 'Error memory allocation dum'
            errnr = 94
            goto 999
         END IF

         do j=1,nanz0
            IF (crtf) THEN
               read(kanal,*) ic(j),ip(j),dum(j)
            ELSE
               read(kanal,*) elec1,elec2,elec3,elec4,dum(j)
               ic(j) = elec1*10000 + elec2
               ip(j) = elec3*10000 + elec4
            END IF
         end do
         close(kanal)

         open(kanal,file=dfm0,status='old')
         read(kanal,*)
         do j=1,nanz0
            read(kanal,*) i,i,dum2(j),idum(j)
         end do
         close(kanal)

         j0 = 0
         i  = 0
 10      i  = i+1
         j = j0
 20      j = j+1
         if (strnr(i).eq.ic(j).and.vnr(i).eq.ip(j).and.
     1        idum(j).eq.1) then
c     nur falls jede Messkonfiguration nur einmal!
c     j0     = j
            d0(i)  = dcmplx(-dlog(dum(j)),0d0)
            fm0(i) = dcmplx(-dlog(dum2(j)),0d0)
         else if (j.lt.nanz0) then
            goto 20
         else
            write(fprun,'(i7,1x,i7,a12)',err=999)
     1           strnr(i),vnr(i),' : discarded'

            nanz = nanz-1
            do j=i,nanz
               strnr(j) = strnr(j+1)
               vnr(j)   = vnr(j+1)
               dat(j)   = dat(j+1)
               wmatd(j) = wmatd(j+1)
               if (lfphai) wmatdp(j)=wmatdp(j+1)
c     nicht notwendig, da Werte eh alle 1
c     wdfak(j) = wdfak(j+1)
            end do
            i = i-1
         end if
         if (i.lt.nanz) goto 10

         open(kanal,file=dm0,status='old')
         read(kanal,*)
         do j=1,elanz
            read(kanal,*) dum3,dum3,dum3
            m0(mnr(j)) = dcmplx(-dlog(1d1)*dum3,0d0)
         end do
         close(kanal)
         DEALLOCATE (dum,dum2,idum,ic,ip)
      end if
c     diff+>
      
      errnr = 0

      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 999  return

 1001 errnr = 2
      return

      end

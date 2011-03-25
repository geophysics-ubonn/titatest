subroutine rall(kanal,delem,delectr,dstrom,drandb,&
!!!$     diff-     1                  dsigma,dvolt,dsens,dstart,lsens,lagain)
!!!$     diff+<
     dsigma,dvolt,dsens,dstart,dd0,dm0,dfm0,lagain)
!!!$     diff+>

!!!$     Unterprogramm zum Einlesen der benoetigten Variablen.

!!!$     Andreas Kemna                                            01-Mar-1995
!!!$     Letzte Aenderung   20-Aug-2007
!!!$     
!!!$.....................................................................

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


!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Kanalnummer
  INTEGER (KIND = 4) ::     kanal

!!!$     Dateinamen
  CHARACTER (80) :: delem,delectr,dstrom,dsigma,dvolt,dsens,dstart, &
!!!$     diff+<
       dd0,dm0,dfm0,drandb
!!!$     diff+>

!!!$     Schalter ob weiterer Datensatz invertiert werden soll
  LOGICAL ::     lagain
  LOGICAL ::     lsto
!!!$     check whether the file format is crtomo konform or not..
  LOGICAL ::    crtf
!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Indexvariable
  INTEGER (KIND = 4) ::     i

!!!$     Pi
  REAL(KIND(0D0))   ::     pi

!!!$     diff+<
  REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: dum,dum2
  REAL(KIND(0D0))                            :: dum3
  INTEGER(KIND = 4),DIMENSION(:),ALLOCATABLE :: idum,ic,ip
  INTEGER (KIND = 4) ::  nanz0,j,j0
!!!$     diff+>

!!!$     ak Inga
  INTEGER (KIND = 4) ::  elec1,elec2,elec3,elec4
!!!$.....................................................................

  pi = dacos(-1d0)

!!!$     'crtomo.cfg' EINLESEN
  fetxt = 'crtomo.cfg'
  errnr = 3
!!!$#################DEFAULTS ########################
!!!$###  switches 
!!!$     ro        lsr    = .false.
!!!$     ro        lpol   = .true.
!!!$     ro        lfphai = .true.
!!!$     ro        lrho0  = .false.
!!!$     akERT2003
!!!$     ak        ld!!!$    = .true.
!!!$     ak        lsr    = .false.
!!!$     ak        lpol   = .false.
!!!$     Sonstiges
!!!$     diff+<
  lsr    = .false.
  lpol   = .false.
  lindiv = .false.
!!!$     ak
!!!$     'dstart'
  lstart = .false.
  dstart = ' '
!!!$     ak        lstart = .true.
!!!$     ak        dstart = '..\..\strasbrg\9610\plane45\mod\rho0.dat'
  ltri   = 0
  lsto = .false.            !default--
!!!$     "Force negative phase" ?
!!!$     sandra        lphi0 = .true.
  lphi0 = .FALSE.
!!!$     ak        lphi0 = .false.
!!!$     "ratio-dataset" ?
  lratio = .false.
!!!$     ak        lratio = .true.
  llamf = .FALSE.
!!!$     final phase improvement setzt phase zueruck auf homogenes modell
  lffhom = .FALSE.
!!!$     Daten Rauschen vom Fehlermodell entkoppeln ?
  lnse2 = .FALSE.
!!!$     Regularisierung mit prior modell?
  lprior = .FALSE.
!!!$######values..
!!!$     FIXED PARAMETER
!!!$     Slash
  slash = '/'
!  CALL CALL get_environment_variable('DELIMITER',slash) ! seems a special C extension 
!!!$     Minimale "L1-ratio" (Grenze der "robust inversion")
  l1min = 1d0
!!!$     ak        l1min = 1.2d0
  nrmsdm = 1d0
!!!$     Art der Ruecktransformation
!!!$     ak        swrtr = 1
!!!$     Minimaler Quotient zweier aufeinanderfolgender Daten-RMS-Werte
!!!$     ak Default
!!!$     ak        mqrms = 1d-2
!!!$     ak ERT 2002-2003 synth
  mqrms = 2d-2
!!!$     ak        mqrms = 2d-2
!!!$     ak Tank
!!!$     ak        mqrms = 2d-2
!!!$     ak MMAJ
!!!$     ak        mqrms = 5d-2
!!!$     CG-Epsilon
  eps = 1d-4
  print*
  print*,eps
  print*
!!!$     Mindest-step-length
  stpmin = 1d-3
!!!$     Minimale stepsize (bdpar)
  bdmin = 0.0d0
!!!$     Regularisierungsparameter
!!!$     ak Default
  nlam   = 30
!!!$     ak Default
  fstart = 0.5d0
  fstop  = 0.9d0
  lamfix = 0.0D0
!!!$     ak MMAJ
!!!$     ak        fstart = 0.2d0
!!!$     ak        fstop  = 0.8d0
!!!$     ak Strasbrg/Werne/Grimberg
!!!$     ak        fstart = 0.5d0
!!!$     ak        fstop  = 0.8d0
  iseedpri = 0; modl_stdn = 0.; iseed = 1
  mswitch = 0
!!!$#########################################################
!!!$     Read in input values..

  fetxt = 'rall -> mswitch'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=98) mswitch

98 fetxt = 'rall -> grid file'
  CALL read_comments(fpcfg)
  READ (fpcfg,'(a)',end=1001,err=999) delem

  fetxt = 'rall -> electrode file'
  CALL read_comments(fpcfg)
  READ (fpcfg,'(a)',end=1001,err=999) delectr

  fetxt = 'rall -> meaurement file'
  CALL read_comments(fpcfg)
  READ (fpcfg,'(a)',end=1001,err=999) dstrom

  fetxt = 'rall -> directory for inversion results'
  CALL read_comments(fpcfg)
  READ (fpcfg,'(a)',end=1001,err=999) ramd
!!$! checks if dir exists and if not, create it
!#if defined (__INTEL_COMPILER)
!!$! check for the intel compiler..
!#define macro_1  INQUIRE ( DIRECTORY=TRIM(ramd),EXIST= crtf)
!#else   
!!$! other compilers go here
!!$! here we may put #elif defined (__GFORTRAN__) as well
!#define macro_1  INQUIRE ( FILE=TRIM(ramd),EXIST= crtf)
!#endif
! ifort uses DIRECTORY for folders, so this is to be used than..
! INQUIRE ( DIRECTORY=TRIM(ramd),EXIST= crtf)

!!$! workaround for compability issues with ifort..
  fetxt = TRIM(ramd)//slash//'tmp.check'
  crtf = .FALSE.
!!$ test if you can open a file in the directory..
  OPEN (fprun,FILE=TRIM(fetxt),STATUS='replace',ERR=97)
  !!$ if you can, the directory exits and you can remove it safely
  CLOSE(fprun,STATUS='delete')
!!$ set this switch to circumvent mkdir
  PRINT*,'writing inversion results into '//TRIM(ramd)
  crtf = .TRUE.
97 IF (.NOT.crtf) THEN
     PRINT*,'Creating inversion directory '//TRIM(ramd)
     CALL SYSTEM ('mkdir '//TRIM(ramd))
  END IF
  fetxt = 'rall -> Difference inversion ?'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) ldiff

  fetxt = 'rall -> Diff. measurements'
  CALL read_comments(fpcfg)
  READ (fpcfg,'(a80)',end=1001,err=999) dd0

  fetxt = 'rall -> Diff. model (prior)'
  CALL read_comments(fpcfg)
  READ (fpcfg,'(a80)',end=1001,err=999) dm0

  fetxt = 'rall -> Diff. model response of prior'
  CALL read_comments(fpcfg)
  READ (fpcfg,'(a80)',end=1001,err=999) dfm0

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
!!!$     diff+>
  fetxt = 'trying noise model seed'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=99) iseedpri,modl_stdn
  !     hier landet man nur, wenn man iseed und modl_stdn angenommen hat
  !      lnse2 = .NOT.lprior       ! kein prior?
  !     Daten Rauschen unabhängig vom Fehlermodell?
  lnsepri = lprior          ! if we have seed and std we assume to add noise to prior
99 fetxt = 'rall -> Gitter nx'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) nx
  fetxt = 'rall -> Gitter nz'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) nz
  fetxt = 'rall -> Anisotropie /x'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) alfx
  fetxt = 'rall -> Anistotropie /y'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) alfz
  fetxt = 'rall -> Maximale Iterationen'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) itmax
!!!$     ak        READ (fpcfg,*,end=1001,err=999) nrmsdm
  fetxt = 'rall -> DC/IP Inversion'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) ldc
!!!$     ak        READ (fpcfg,*,end=1001,err=999) lsr
  fetxt = 'rall -> Robuste Inversion'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) lrobust
  IF (lrobust) PRINT*,'## Robust inversion ##'
!!!$     ak        READ (fpcfg,*,end=1001,err=999) lpol
  fetxt = 'rall -> Finale Phasen Inversion'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) lfphai
!!!$     ak        READ (fpcfg,*,end=1001,err=999) lindiv
  fetxt = 'rall -> Relativer Fehler Widerstand'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) stabw0
  fetxt = 'rall -> Absoluter Fehler Widerstand'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) stabm0
  fetxt = 'rall -> Phasenfehlerparameter A1'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) stabpA1
  fetxt = 'rall -> Phasenfehlerparameter B'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) stabpB
  fetxt = 'rall -> Relative Fehler Phasen'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) stabpA2
  fetxt = 'rall -> Absoluter Fehler Phasen (mRad)'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) stabp0
  fetxt = 'rall -> Homogenes Startmodell?'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) lrho0
  fetxt = 'rall -> rho_0'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) bet0
  fetxt = 'rall -> phase_0'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) pha0
  fetxt = 'rall -> Noch eine Inversion'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) lagain
  fetxt = 'rall -> 2D oder 2.5D ?'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) swrtr
  fetxt = 'rall -> weitere Quelle?'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) lsink
  fetxt = 'rall -> Nummer der Quelle'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) nsink
  fetxt = 'rall -> Randbedingungen ?'
  CALL read_comments(fpcfg)
  READ (fpcfg,*,end=1001,err=999) lrandb2
  fetxt = 'rall -> Datei mit Randwerten'
  CALL read_comments(fpcfg)
  READ (fpcfg,'(a80)',end=1001,err=999) drandb
  fetxt = 'triangularization switch'
  CALL read_comments(fpcfg)
  READ (fpcfg,'(I2)',end=100,err=100) ltri

  IF (BTEST(ltri,5)) THEN
     llamf = .TRUE.
     fetxt = 'rall -> fixed lam value'
     CALL read_comments(fpcfg)
     READ (fpcfg,*,end=104,err=104) lamfix
     GOTO 105
104  lamfix = 1.0           ! default value for MGS
     BACKSPACE(fpcfg)

105  PRINT*,'Fixing Lambda =', lamfix
     ltri = ltri - 2**5
  END IF

  IF (ltri > 15) THEN ! exception for wrong ltri switch
     PRINT*,'WARNING, fix lambda switch has changed'
     PRINT*,'check ltri value (>15):',ltri
     STOP
  END IF

  lsto = (ltri == 15)

  GOTO 101

100 BACKSPACE (fpcfg)


101 IF (lsto) PRINT*,'Stochastische Regularisierung'

  IF (ltri > 4 .AND. ltri < 15) THEN
     fetxt = 'rall -> beta value'
     CALL read_comments(fpcfg)
     READ (fpcfg,*,end=102,err=102) betamgs
     GOTO 103
102  betamgs = 0.1          ! default value for MGS
     BACKSPACE (fpcfg)

103  PRINT*,'Regularisation with support stabilizer beta =',betamgs
  END IF

  IF (itmax == 0) PRINT*,' ####### Only precalcs, itmax==0 ###########'

!!!$     check if the final phase should start with homogenous model      
  lffhom = (stabp0 < 0)
  IF (lffhom) stabp0 = -stabp0

!!!$     check if there is crt.noisemod containig noise info
  fetxt = 'crt.noisemod'
  INQUIRE(FILE=TRIM(fetxt),EXIST=lnse2)

  lnse = ( stabw0 < 0 )     ! couple error and noise model
  IF ( lnse ) THEN
     stabw0 = -stabw0
     IF (lnse2) print*,'overriding seperate noise model'
     lnse2 = .FALSE.        ! overrides the lnse2 switch
!!!$     copy error model into noise model
     nstabw0 = stabw0
     nstabm0 = stabm0
     nstabpA1 = stabpA1
     nstabpA2 = stabpA2
     nstabp0 = stabp0

     fetxt = 'rall -> seed'
     CALL read_comments(fpcfg)
     READ (fpcfg,*,end=106,err=106) iseed
     GOTO 107
106  iseed = 1              ! default value for PRS
     BACKSPACE(fpcfg)
     WRITE (*,'(a)')' Rauschen Gekoppelt an Fehlermodell '
  END IF

107 IF (lnse2) THEN

     fetxt = 'get noise model from crt.noisemod'
     CALL get_noisemodel(iseed,nstabw0,nstabm0,nstabpA1,&
          nstabpB,nstabpA2,nstabp0,errnr)

     IF (errnr /= 0) GOTO 999

     WRITE (*,'(a,I7)',ADVANCE='no')&
          'Entkoppeltes Daten Rauschen:: seed:',iseed

     lnse = .TRUE.          ! add noise

  END IF

  IF (lnse) THEN 
     fetxt = 'write out noise model'
     CALL write_noisemodel(iseed,nstabw0,nstabm0,&
          nstabpA1,nstabpB,nstabpA2,nstabp0,errnr)
     IF (errnr /= 0) GOTO 999
  ELSE
     PRINT*,'No Data noise!!'
  END IF

  IF ((nx<=0.OR.nz<=0).AND.ltri==0) ltri=1 ! at least L1-smoothness

!!!$     Ggf. Fehlermeldungen
  if (ltri==0.AND.(nx.lt.2.or.nz.lt.2)) then
     fetxt = ' '
     errnr = 89
     goto 999
!!!$  else if (alfx.le.0d0.or.alfz.le.0d0) then
!!!$  fetxt = ' '
!!!$  errnr = 96
!!!$  goto 999
  else if (itmax<0.or.itmax.ge.100) then
     fetxt = ' '
     errnr = 61
     goto 999
  else if (nrmsdm.lt.1d-12) then
     fetxt = ' '
     errnr = 62
     goto 999
!!!$     else if (nlam.lt.0) then
!!!$     fetxt = ' '
!!!$     errnr = 83
!!!$     goto 999
!!!$     else if (fstart.gt.1d0.or.fstop.gt.1d0.or.
!!!$     1           fstart.le.0d0.or.fstop.le.0d0.or.
!!!$     1           fstart.gt.fstop) then
!!!$     fetxt = ' '
!!!$     errnr = 98
!!!$     goto 999
  else if (stabw0.le.0d0.or.stabm0.lt.0d0) then
     fetxt = ' '
     errnr = 104
     goto 999
  else if (.not.ldc.and.lfphai.and.&
    ((stabp0.lt.0d0.or.stabpA2.lt.0d0).OR. &
    ((stabp0 == 0d0).and.(stabpA2 == 0d0)))) then
     fetxt = ' '
     errnr = 105
     goto 999
  else if (lrho0.and.(bet0.le.0d0.or.&
       (.not.ldc.and.dabs(pha0).gt.1d3*pi))) then
     fetxt = ' '
     errnr = 91
     goto 999
!!!$     else if (mqrms.lt.0d0.or.mqrms.ge.1d0) then
!!!$     fetxt = ' '
!!!$     errnr = 64
!!!$     goto 999
!!!$     else if (lrobust.and.l1min.lt.1d0) then
!!!$     fetxt = ' '
!!!$     errnr = 90
!!!$     goto 999
  end if

  lelerr = .NOT.lfphai.AND..NOT.ldc ! complex inversion only

!!!$     (mswitch) Mega switch testing..
  lsens = BTEST(mswitch,0)  ! +1 ueberdeckung schreiben
  lcov1 = BTEST(mswitch,1)  ! +2 posterior modell covariance matrix 1
  lres  = BTEST(mswitch,2)  ! +4 rsolution matrix berechnen
  lcov2 = BTEST(mswitch,3)  ! +8 posterior modell covariance matrix 2

  lgauss = BTEST (mswitch,4) ! +16 solve ols with Gauss elemination

  lelerr = BTEST (mswitch,5).OR.lelerr ! +32 overwrites previous lelerr

  IF (lelerr) PRINT*,'## using complex error ellipses ##'

  lphi0 = BTEST (mswitch,7) ! +128 forcing negative phase

  lsytop = .NOT.BTEST (mswitch,8) ! +256 disables sy top check of 
  !     no flow boundary electrodes for enhanced beta calculation (bsytop). 
  !     This is useful for including topographical effects and should be used

  lvario = BTEST (mswitch,9) ! +512 calculate variogram

  lverb = BTEST (mswitch,10) ! +1024 Verbose output CG, daten, bnachbar..

  IF (lverb) WRITE(*,'(/a/)')' #  ## VERBOSE ## #'

  lres = (lres.or.lcov2)    ! compute mcm2 on top of resolution
  lcov1 = (lres.or.lcov1)   ! compute resolution by taking mcm1
!!!$     
  lsens = .TRUE.            ! default immer coverages schreiben..
!!!$     
  if (lratio) then
     lrho0  = .true.
     lstart = .false.
     lphi0  = .false.
     lpol   = .false.
  end if
!!!$     diff-        if (lstart) lrho0=.false.

  if (lstart.or.ldiff) lrho0=.false.
!!!$     ak
  if (ldiff) then
     ldc  = .true.
     lpol = .false.
  end if
!!!$     diff+>
!!!$     ak        if (ldc.or.stabp0.ge.stabw0) lfphai=.false.
  if (ldc) lfphai=.false.
!!!$     Dateien
  lnramd = index(ramd,' ')-1
  dsigma = ramd(1:lnramd)//slash(1:1)//'rho.dat'
  dvolt  = ramd(1:lnramd)//slash(1:1)//'volt.dat'
  dsens  = ramd(1:lnramd)//slash(1:1)//'coverage.mag'

!!!$     Elementeinteilung einlesen
  IF (lverb) WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//'reading grid'
  call relem(kanal,delem)
  if (errnr.ne.0) goto 999

  IF (ltri/=0) THEN
     manz = elanz           ! wichtig an dieser stelle..
     CALL bnachbar          ! blegt nachbar
     CALL besp_elem
     lvario = lvario.OR.lsto
  ELSE
!!!$     Modelleinteilung gemaess Elementeinteilung belegen
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
  lvario = lvario.OR. &       ! if already set or
       (itmax == 0).AND.(lstart.OR.lprior) ! analyse any prior

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

!!!$     Maximale Anzahl an CG-steps setzen
!!!$     ak        ncgmax = manz
!!!$  ncgmax = manz / 2         ! useful for small scale model variations
  ncgmax = manz         ! useful for small scale model variations
  !     for normal smooth and damping we usually need fewer CG iterations;
  !     because the model variations are of bigger scale size
  IF (ltri < 5) ncgmax = ncgmax / 10
!!!$     Elektrodenverteilung und Daten einlesen sowie Wellenzahlwerte
!!!$     bestimmen

  call relectr(kanal,delectr)
  if (errnr.ne.0) goto 999

  IF (lsink) THEN
     IF (nsink > sanz) THEN
        PRINT*,'Sink node > grid nodes'
        errnr = 3
        GOTO 999
     END IF
     WRITE(*,'(/A,I5,2F12.3/)')'Fictious sink @ node ',&
          nsink,sx(snr(nsink)),sy(snr(nsink))
!!!$         WRITE(fpinv,'(A,I5,2F12.3)')'Fictious sink @ node ',
!!!$     1        nsink,sx(snr(nsink)),sy(snr(nsink))
  END IF

  call rdati (kanal,dstrom)

  if (errnr.ne.0) goto 999

  if (swrtr.eq.0) then
     lsr    = .false.
     kwnanz = 1
     ALLOCATE (kwn(kwnanz),kwnwi(kwnanz),stat=errnr)
     IF (errnr /= 0) THEN
        fetxt = 'Error memory allocation kwn'
        errnr = 94
        GOTO 999
     END IF
     kwn = 0d0; kwnwi = 0D0
     do i=1,typanz
        IF (typ(i) == 11) THEN
           fetxt = 'in 2D keine gemischten RB'
           PRINT*,TRIM(fetxt)
           errnr = 110
           GOTO 999
        END IF
     END DO
     IF (.NOT.lrandb2) THEN
        WRITE (*,'(//a,t33,a)')'2D without Dirichlet nodes','setting k=1e-6'
        kwn = 1d-6  ! make sure A is still pos definite
     END IF
  else
     call rwaven()
     if (errnr.ne.0) goto 999
  end if

!!!$     read boundary values
  if (lrandb2) then
     call rrandb(kanal,drandb)
     if (errnr.ne.0) goto 999
  end if
!!!$     diff+<
  if (ldiff) then
     ALLOCATE (d0(nanz),fm0(nanz),stat=errnr)
     IF (errnr /= 0) THEN
        fetxt = 'Error memory allocation diff data '
        errnr = 94
        goto 999
     END IF
     open(kanal,file=TRIM(dd0),status='old')
     read(kanal,*) nanz0
     read(kanal,*,err=999) elec1
     BACKSPACE(kanal)

     elec3=elec1-10000      ! are we still positive?
     crtf=(elec3 > 0)       ! crtomo konform?

     ALLOCATE (dum(nanz0),dum2(nanz0),idum(nanz0),&
          ic(nanz0),ip(nanz0),stat=errnr)
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

     open(kanal,file=TRIM(dfm0),status='old')
     read(kanal,*)
     do j=1,nanz0
        read(kanal,*) i,i,dum2(j),idum(j)
     end do
     close(kanal)

     j0 = 0
     i  = 0
10   i  = i+1
     j = j0
20   j = j+1
     if (strnr(i).eq.ic(j).and.vnr(i).eq.ip(j).and.idum(j).eq.1) then
!!!$     nur falls jede Messkonfiguration nur einmal!
!!!$     j0     = j
        d0(i)  = dcmplx(-dlog(dum(j)),0d0)
        fm0(i) = dcmplx(-dlog(dum2(j)),0d0)
     else if (j.lt.nanz0) then
        goto 20
     else
        write(fprun,'(i7,1x,i7,a12)',err=999)strnr(i),vnr(i),' : discarded'

        nanz = nanz-1
        do j=i,nanz
           strnr(j) = strnr(j+1)
           vnr(j)   = vnr(j+1)
           dat(j)   = dat(j+1)
           wmatd(j) = wmatd(j+1)
           if (lfphai) wmatdp(j)=wmatdp(j+1)
!!!$     nicht notwendig, da Werte eh alle 1
!!!$     wdfak(j) = wdfak(j+1)
        end do
        i = i-1
     end if
     if (i.lt.nanz) goto 10

     open(kanal,file=TRIM(dm0),status='old')
     read(kanal,*)
     do j=1,elanz
        read(kanal,*) dum3,dum3,dum3
        m0(mnr(j)) = dcmplx(-dlog(1d1)*dum3,0d0)
     end do
     close(kanal)
     DEALLOCATE (dum,dum2,idum,ic,ip)
  end if ! ldiff
!!!$     diff+>

  errnr = 0

  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

999 return

1001 errnr = 2
  return

end subroutine rall


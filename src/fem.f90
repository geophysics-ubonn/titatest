program fem
  ! Complex Resistivity 2.5D Modelling Main Program
  !
  ! Andreas Kemna                                            15-Dec-1994
  ! Letzte Aenderung   16-Sep-2014
  !....................................................................

  use alloci
  use tic_toc
  use femmod
  use datmod
  use sigmamod
  use electrmod
  use modelmod
  use elemmod
  use wavenmod
  use randbmod
  use errmod
  use konvmod
  use omp_lib
  use ompmod
  use get_ver

  implicit none

  !....................................................................

  ! channel number
  integer               kanal

  ! file names
  character (len = 80)  delem,delectr,dsigma,dstrom,&
       dkpot,dpot,dvolt,dsens,drandb

  ! switch write Fourier-space potentials
  logical               lkpot

  ! switch write normal-space potentials
  logical               lpot

  ! switch write voltages/measurements
  logical               lvolt

  ! switch run program again
  logical               lagain

  ! switch modelling with K factor (default .false.)
  logical               wkfak

  ! switch analytical modelling (default .false.)
  logical               lana

  ! index variables
  integer               j, k, l, c1

  ! error message variable
  character(len = 256)  ftext
  integer               getpid, pid 

  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  kanal = 11
  pid = getpid()
  fetxt = 'crmod.pid'
  print*,'CRMOD complex resistivity modelling'
  print*
  print*,'program info'
  print*,'crmod pid  ',pid
  open (kanal, FILE = trim(fetxt), STATUS = 'replace', err = 999)
  write (kanal,*) pid
  close (kanal)

  call get_git_ver( version )

  call tic( c1 )
  ! open 'crmod.cfg' 
  fetxt = 'crmod.cfg'
  errnr = 1
  open (12, file = trim( fetxt ), status = 'old', err = 999)

  ! set general parameters
  ! singularity removal
  lsr = .false.

  ! other
  lkpot = .false.
  wkfak = .false.
  lana = .false.
  lprior = .false. ! kann prior modell verrauschen
  ldc = .false.
  lnsepri = .false.

  ! control output
5 write (*,*)
  write (*,'(a20)') ' reading input files'

  ! read 'crmod.cfg'
  fetxt = 'crmod.cfg'
  errnr = 3
  mswitch = 0
  read (12, *,       end = 1001, err = 999)
  read (12, '(a80)', end = 1001, err = 999) delem
  read (12, '(a80)', end = 1001, err = 999) delectr
  read (12, '(a80)', end = 1001, err = 999) dsigma
  read (12, '(a80)', end = 1001, err = 999) dstrom
  read (12, *,       end = 1001, err = 999) lpot
  read (12, '(a80)', end = 1001, err = 999) dpot
  read (12, *,       end = 1001, err = 999) lvolt
  read (12, '(a80)', end = 1001, err = 999) dvolt
  read (12, *,       end = 1001, err = 999) lsens
  read (12, '(a80)', end = 1001, err = 999) dsens
  read (12, *,       end = 1001, err = 999) lagain
  read (12, *,       end = 1001, err = 999) swrtr
  read (12, *,       end = 1001, err = 999) lsink
  read (12, *,       end = 1001, err = 999) nsink
  read (12, *,       end = 1001, err = 999) lrandb2
  read (12, '(a80)', end = 1001, err = 999) drandb
  fetxt = 'reading mswitch'
  read (12, '(I4)',  end = 101,  err = 100) mswitch

  goto 101

100 backspace (12)

101 lana = btest(mswitch,0)
  wkfak = btest(mswitch,1)
  lsr = btest(mswitch,2)
  lverb = btest(mswitch,10)

  if (lana) print*,'analytical modelling'
  if (wkfak) print*,'modelling with K factor'
  if (lsr) then
     print*,'modelling with sigularity removal'
     lana = .false. ! analytical for singularity
     ! is indeed controlled by lsr and not lana..
     ! lana is only true for analytical solution only
  end if

  if (lverb) print*,'verbose output'

  ! read grid: element file
  call relem(kanal,delem)
  if (errnr.ne.0) goto 999

  ! read grid: electrode positions
  call relectr(kanal,delectr)
  if (errnr.ne.0) goto 999

  ! read grid: electrode configurations
  call rdatm(kanal,dstrom)
  if (errnr.ne.0) goto 999

  ! 2D: swrtr.eq.0; 2.5D: swrtr.eq.1
  if (swrtr.eq.0) then
     ! no analytical modelling and no singularity removal in 2D
     lana = .false.
     lsr = .false.
     ! only one wavenumber with value 0. this transforms the 2.5D Helmholtz
     ! equation into a 2D Poisson equation
     kwnanz = 1
     allocate (kwn(kwnanz),stat=errnr)
     if (errnr /= 0) then
        fetxt = 'error memory allocation kwn'
        errnr = 94
        goto 999
     end if
     kwn(1) = 0d0
  else
     ! compute wave numbers for 2.5D modelling
     call rwaven()
     if (errnr.ne.0) goto 999
  end if

  allocate (sigma(elanz),stat=errnr)
  if (errnr /= 0) then
     fetxt = 'error memory allocation fem sigma'
     errnr = 94
     goto 999
  end if

  manz =elanz

  ! read grid: cell parameters
  call rsigma(kanal,dsigma)
  if (errnr.ne.0) goto 999

  ! read fixed boundary values
  if (lrandb2) then
     call rrandb(kanal,drandb)
     if (errnr.ne.0) goto 999
  end if

  ! get grid-averaged conductivity
  call refsig()

  ! print fictious sink node and position
  if (lsink) write(6,'(A,I6,A,F10.2,A,F10.2,A)')&
       ' fictious sink at node',nsink,', x=',sx(snr(nsink)),&
       'm, y=',sy(snr(nsink)),'m'
  ! pre-assemble stiffness matrices from area elements
  ! in case of mixed boundaries (ntyp.eq.8), compute beta value 
  call precal()
  if (errnr.ne.0) goto 999

  if (.not.lbeta) lsr=.false.

  if (.not.allocated(hpot)) allocate(hpot(sanz,eanz))
  if (.not.allocated(kpot)) allocate(kpot(sanz,eanz,kwnanz))
  if (.not.allocated(b)) allocate(b(sanz,1))
  if (.not.allocated(a)) allocate(a(2*mb+1,sanz))
  !  if (.not.allocated(a_mat_band)) allocate(a_mat_band(2*mb+1,sanz))
  if (.not.allocated(x)) allocate(x(sanz,1))
  !  if (.not.allocated(a_mat_band_elec)) allocate(a_mat_band_elec(mb+1,sanz))
  if (.not.allocated(pot))  allocate(pot(sanz))
  if (.not.allocated(pota))  allocate(pota(sanz))
  if (.not.allocated(fak))  allocate(fak(sanz))
  print*
  print*,'finite element modelling'
  ! compute potentials
  do k=1,kwnanz
     if (swrtr.eq.0) then
        write(*,'(a)')' computing potentials'
     else
        write (*,'(a,I4,a,I4)',ADVANCE='no')&
             achar(13)//' computing potentials : wavenumber',k,' of',kwnanz
     end if
     do l=1,eanz
        if (lsr.or.lbeta.or.l.eq.1) then
!   Ggf. Potentialwerte fuer homogenen Fall analytisch berechnen
           if (lsr) call potana(l,k,pota)

!   Kompilation des Gleichungssystems (fuer Einheitsstrom !)
           call kompab(l,k)
           if (errnr.ne.0) goto 999

!   Ggf. Randbedingung beruecksichtigen
           if (lrandb) call randb()
           if (lrandb2) call randb2()

        else
!   Stromvektor modifizieren
           call kompb(l)
        end if

!   Gleichungssystem loesen
       IF (.NOT.lana) call solve_zgbsvx(a,x,-b)
!   Potentialwerte zurueckskalieren und umspeichern sowie ggf.
!   analytische Loesung addieren
        do j=1,sanz
           if (lana) then
              kpot(j,l,k) = pota(j)
           else
              kpot(j,l,k) = x(j,1) 
              if (lsr) kpot(j,l,k) = kpot(j,l,k) + pota(j)
           end if
           if (swrtr.eq.0) hpot(j,l) = kpot(j,l,k)
        end do
     end do
  end do
  print*,''
  ! inverse Fourier transform
  if (swrtr.eq.1) call rtrafo()
  print*
  print*,'writing output files'
  ! write potentials in wavenumber space (if(lkpot))
  if (lkpot) then
     call wkpot(kanal,dkpot)
     if (errnr.ne.0) goto 999
  end if

  ! compute potentials and write to file (if(lpot))
  if (lpot) then
     call bpot(kanal,dpot)
     if (errnr.ne.0) goto 999
  end if

  ! compute voltages and write to file (if(lvolt))
  if (lvolt) then
     ! compute apparent resistivities (if(wkfak.and.lbeta)). Only possible for flat 
     ! surfaces (configuration factor has to exist)
     if (wkfak.and.lbeta) then
        call bkfak()
        if (errnr.ne.0) goto 999
     else
        kfak=1.0 
     end if

     call bvolt()
     if (errnr.ne.0) goto 999

     call wdatm(kanal,dvolt)
     if (errnr.ne.0) goto 999
  end if

  ! compute sensitivities of all recordings (if(lsens))
  if (lsens) then
     if (manz.ne.elanz) then
        fetxt = 'manz /= elanz .. is not implemented yet'
        errnr = 50
        goto 999
     end if
     allocate (mnr(elanz))
     if (errnr /= 0) then
        fetxt = 'Error memory allocation mnr failed'
        errnr = 94
        goto 999
     end if
     ! set parameter to grid cell relation
     do j=1,elanz
        mnr(j) = j
     end do

     allocate(sens(nanz,manz),stat=errnr)
     if (errnr.ne.0) then
        errnr = 97 
        goto 999
     end if
     ! compute and write sensitivities 
     call bsens()
     call wsens(kanal,dsens)
     if (errnr.ne.0) goto 999
     deallocate(sens)
  end if

  deallocate(kpot,pot,pota,fak,a)

  ! run program again
  if (lagain) goto 5

  ! close 'crmod.cfg'
  close(12)

  ! control output
  fetxt = 'solution time'
  call TOC(c1,fetxt)

  write(*,*)
  if (allocated (snr))   deallocate (snr, sx, sy)
  if (allocated (typ))   deallocate (typ, nelanz, selanz)
  if (allocated (nrel))  deallocate (nrel, rnr)
  if (allocated (strnr)) deallocate (strnr, strom, volt, sigmaa, kfak, vnr)
  if (allocated (sigma)) deallocate (sigma)
  if (allocated (enr))   deallocate (enr)
  if (allocated (mnr))   deallocate (mnr)
  if (allocated (kwn))   deallocate (kwn)
  if (allocated (kwnwi)) deallocate (kwnwi)
  if (allocated (rwddc)) deallocate (rwddc) 
  if (allocated (rwndc)) deallocate (rwndc) 
  if (allocated (rwd))   deallocate (rwd) 
  if (allocated (rwn))   deallocate (rwn) 
  if (allocated (rwdnr)) deallocate (rwdnr) 
  write(*,'(a)',advance='no') ' '
  stop 'modelling completed'

  !....................................................................

  ! write error messages
999 open(9,file = 'error.dat',status = 'replace')
  errflag = 2
  call get_error(ftext,errnr,errflag,fetxt)
  write(9,'(a80,i3,i1)') fetxt,errnr,errflag
  write(9,*)ftext
  close(9)
  stop 'error -1'

1001 open(9,file = 'error.dat',status = 'replace')
  errnr   = 2
  errflag = 2
  call get_error(ftext,errnr,errflag,fetxt)
  write(9,'(a80,i3,i1)') fetxt,errnr,errflag
  write(9,*)ftext
  close(9)
  stop 'error -2'

end program fem

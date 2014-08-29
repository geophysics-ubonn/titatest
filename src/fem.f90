program fem

!   Hauptprogramm zum Complex-Resistivity-2.5D-Modelling.

!   Belegte Kanaele:  9 - error.dat
!   11 - in-/output
!   12 - crmod.cfg

!   Andreas Kemna                                            15-Dec-1994
!   Letzte Aenderung   02-May-2008

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

!   Kanalnummer
  integer             kanal


!   Dateinamen
  character (len=80)  delem,delectr,dsigma,dstrom,&
       dkpot,dpot,dvolt,dsens,drandb

!   Schalter ob transformierte Potentialwerte ausgegeben werden sollen
  logical             lkpot

!   Schalter ob Potentialwerte ausgegeben werden sollen
  logical             lpot

!   Schalter ob Spannungswerte ausgegeben werden sollen
  logical             lvolt

!   Schalter ob weiterer Datensatz modelliert werden soll
  logical             lagain

!   Schalter ob mit K-Faktor modelliert werden soll (default ohne)
  logical             wkfak

!   Schalter ob nur analytisch modelliert werden soll (default ohne)
  logical             lana

!   Indexvariablen
  integer             j,k,l,c1,info

! counting wavenumbers
  integer :: count

  character(len=256)    ftext
  integer  getpid,pid,maxthreads,mythreads

!  ! potential solution vector
!  complex(kind(0d0)),dimension(:,:),allocatable :: x


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  kanal = 11
  pid = getpid()
  fetxt = 'crmod.pid'
  print*,'CRMod Process_ID ::',pid
  open (kanal,FILE=trim(fetxt),STATUS='replace',err=999)
  write (kanal,*)pid
  close (kanal)
  !  maxthreads = OMP_GET_MAX_THREADS()
  !  WRITE(6,"(a, i3)") " OpenMP max threads: ", maxthreads

  call get_git_ver(version)

  call tic(c1)
!   'crmod.cfg' oeffnen
  fetxt = 'crmod.cfg'
  errnr = 1
  open(12,file=trim(fetxt),status='old',err=999)

!   Allgemeine Parameter setzen


!   "singularity removal" ?
!   ak        lsr = .true.
  lsr = .false.

!   Art der Ruecktransformation
!   ak        swrtr = 1
!   ak        swrtr = 0

!   Sonstiges
  lkpot = .false.
  wkfak = .false.
  lana = .false.
  lprior = .false. ! kann prior modell verrauschen
  ldc = .false.
  lnsepri = .false.
!   ak        lkpot = .true.
!   ak        dkpot = '..\tmp\kpot.ctr'

!   Kontrollausgabe
5 write(*,*)
  write(*,'(a20)') ' Reading Input-Files'

!   'crmod.cfg' einlesen
  fetxt = 'crmod.cfg'
  errnr = 3
  mswitch = 0
  read(12,*,end=1001,err=999)
  read(12,'(a80)',end=1001,err=999) delem
  read(12,'(a80)',end=1001,err=999) delectr
  read(12,'(a80)',end=1001,err=999) dsigma
  read(12,'(a80)',end=1001,err=999) dstrom
  read(12,*,end=1001,err=999) lpot
  read(12,'(a80)',end=1001,err=999) dpot
  read(12,*,end=1001,err=999) lvolt
  read(12,'(a80)',end=1001,err=999) dvolt
  read(12,*,end=1001,err=999) lsens
  read(12,'(a80)',end=1001,err=999) dsens
  read(12,*,end=1001,err=999) lagain
  read(12,*,end=1001,err=999) swrtr
  read(12,*,end=1001,err=999) lsink
  read(12,*,end=1001,err=999) nsink
  read(12,*,end=1001,err=999) lrandb2
  read(12,'(a80)',end=1001,err=999) drandb
  fetxt = 'reading mswitch'
  read(12,'(I4)',end=101,err=100) mswitch

  goto 101

100 backspace (12)

101 lana = btest(mswitch,0)
  wkfak = btest(mswitch,1)
  lsr = btest(mswitch,2)
  lverb = btest(mswitch,10)

  if (lana) print*,'ANALYTICAL SOLUTION'
  if (wkfak) print*,'WITH K-FAKTOR'
  if (lsr) then
     print*,'WITH SINGULARITY REMOVAL'
     lana = .false. ! analytical for singularity
     ! is indeed controlled by lsr and not lana..
     ! lana is only true for analytical solution only
  end if

  if (lverb) print*,'VERBOSE OUTPUT'

!   Alles einlesen
  call relem(kanal,delem)
  if (errnr.ne.0) goto 999

  call relectr(kanal,delectr)
  if (errnr.ne.0) goto 999

  call rdatm(kanal,dstrom)
  if (errnr.ne.0) goto 999


  if (swrtr.eq.0) then
     lana = .false.
     lsr = .false.
     kwnanz = 1
     allocate (kwn(kwnanz),stat=errnr)
     if (errnr /= 0) then
        fetxt = 'Error memory allocation kwn'
        errnr = 94
        goto 999
     end if
     kwn(1) = 0d0
  else
     call rwaven()
     if (errnr.ne.0) goto 999
  end if
  print*,''
  print*
!   Startmodell belegen
  allocate (sigma(elanz),stat=errnr)
  if (errnr /= 0) then
     fetxt = 'Error memory allocation fem sigma'
     errnr = 94
     goto 999
  end if

  manz =elanz

  call rsigma(kanal,dsigma)
  if (errnr.ne.0) goto 999

  if (lrandb2) then
     call rrandb(kanal,drandb)
     if (errnr.ne.0) goto 999
  end if

!   Ggf. Referenzleitfaehigkeit bestimmen
  call refsig()

  if (lsink) write(6,'(/A,I5,2F12.3/)')&
       'Fictious sink @ node ',nsink,sx(snr(nsink)),sy(snr(nsink))
!   Element- und Randelementbeitraege sowie Konfigurationsfaktoren zur
!   Berechnung der gemischten Randbedingung bestimmen
  call precal()
  if (errnr.ne.0) goto 999

  if (.not.lbeta) lsr=.false.

  if (.not.allocated(hpot)) allocate(hpot(sanz,eanz))
  if (.not.allocated(kpot)) allocate(kpot(sanz,eanz,kwnanz))
  if (.not.allocated(b)) allocate(b(sanz,1))
  if (.not.allocated(a_mat_band)) allocate(a_mat_band(2*mb+1,sanz))
if (.not.allocated(x)) allocate(x(sanz,1))
if (.not.allocated(a_mat_band_elec)) allocate(a_mat_band_elec(mb+1,sanz))
if (.not.allocated(pot))  allocate(pot(sanz))
if (.not.allocated(pota))  allocate(pota(sanz))
if (.not.allocated(fak))  allocate(fak(sanz))
  count = 0
! Compute potentials
  do k=1,kwnanz
     count = count + 1
     if (swrtr.eq.0) then
        write(*,'(a)')' Computing Potentials'
     else
        write (*,'(a,t45,I4,t100,a)',ADVANCE='no')&
             achar(13)//' Computing Potentials : Wavenumber ',count
     end if
     call pre_comp_ab(k,a_mat_band)
     do l=1,eanz
        b = cmplx(0.)
        a_mat_band_elec = a_mat_band
        b(enr(l),1) = cmplx(1.)
        if (lsink) b(nsink,1) = cmplx(-1.)
        call comp_ab(k,a_mat_band_elec,l)
! Incorporate special boundary conditions
        !           if (lrandb) call randb(a,b)
        !           if (lrandb2) call randb2(a,b)

        ! General Band matrix, expert solver
        call solve_zgbsvx(a_mat_band_elec,x,b)

! Save potentials for each wavenumber (if necess. = if(swrtr)) for Fourier 
! back-transform.
! Add analytical potentials (if available)
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
  print*,'done, now processing'
! Fourier back-transform
  if (swrtr.eq.1) call rtrafo()

! Write potentials in wavenumber space (if(lkpot))
  if (lkpot) then
     call wkpot(kanal,dkpot)
     if (errnr.ne.0) goto 999
  end if

! Compute potentials and write to file (if(lpot))
  if (lpot) then
     call bpot(kanal,dpot)
     if (errnr.ne.0) goto 999
  end if

! Compute voltages and write to file (if(lvolt))
  if (lvolt) then
! Compute apparent resistivities (if(wkfak.and.lbeta)). Only possible for flat 
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

! Compute sensitivities of all recordings (if(lsens))
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
!   Modelleinteilung gemaess Elementeinteilung belegen

     do j=1,elanz
        mnr(j) = j
     end do

!   'sens' zuweisen
     allocate(sens(nanz,manz),stat=errnr)
     if (errnr.ne.0) then
        errnr = 97 
        goto 999
     end if

     call bsens()

     call wsens(kanal,dsens)
     if (errnr.ne.0) goto 999

!   'sens' freigeben
     deallocate(sens)
  end if

!   'kpot' freigeben
  deallocate(kpot,pot,pota,fak,a_mat_band,a_mat_band_elec)

!   Ggf. weiteren Datensatz modellieren
  if (lagain) goto 5

!   'crmod.cfg' schliessen
  close(12)

!   Kontrollausgabe

  fetxt = 'solution time'
  call TOC(c1,fetxt)

  write(*,*)
  write(*,'(a)',ADVANCE='no')' Modelling completed'


  if (allocated (snr)) deallocate (snr,sx,sy)
  if (allocated (typ)) deallocate (typ,nelanz,selanz)
  if (allocated (nrel)) deallocate (nrel,rnr)

  if (allocated (strnr)) deallocate (strnr,strom,volt,sigmaa,&
       kfak,vnr)
  if (allocated (sigma)) deallocate (sigma)
  if (allocated (enr)) deallocate (enr)
  if (allocated (mnr)) deallocate (mnr)
  if (allocated (kwn)) deallocate (kwn)
  if (allocated (kwnwi)) deallocate (kwnwi)

  if (allocated (rwddc)) deallocate (rwddc) 
  if (allocated (rwndc)) deallocate (rwndc) 
  if (allocated (rwd)) deallocate (rwd) 
  if (allocated (rwn)) deallocate (rwn) 
  if (allocated (rwdnr)) deallocate (rwdnr) 

  stop '0'

!....................................................................

!   (Fehler-) Meldung schreiben
999 open(9,file='error.dat',status='replace')
  errflag = 2
  call get_error(ftext,errnr,errflag,fetxt)
  write(9,'(a80,i3,i1)') fetxt,errnr,errflag
  write(9,*)ftext
  close(9)
  stop '-1'

1001 open(9,file='error.dat',status='replace')
  errnr   = 2
  errflag = 2
  call get_error(ftext,errnr,errflag,fetxt)
  write(9,'(a80,i3,i1)') fetxt,errnr,errflag
  write(9,*)ftext
  close(9)
  stop '-2'

end program fem

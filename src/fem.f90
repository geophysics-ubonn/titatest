program fem

!!!$   Hauptprogramm zum Complex-Resistivity-2.5D-Modelling.

!!!$   Belegte Kanaele:  9 - error.dat
!!!$   11 - in-/output
!!!$   12 - crmod.cfg

!!!$   Andreas Kemna                                            15-Dec-1994
!!!$   Letzte Aenderung   02-May-2008

!!!$....................................................................

    USE alloci
    USE tic_toc
    USE femmod
    USE datmod
    USE sigmamod
    USE electrmod
    USE modelmod
    USE elemmod
    USE wavenmod
    USE randbmod
    USE errmod
    USE konvmod
    USE omp_lib
    USE ompmod
    USE get_ver

    IMPLICIT none


    !!!$....................................................................

    !!!$   Kanalnummer
    integer             kanal


    !!!$   Dateinamen
    character (len=80)  delem,delectr,dsigma,dstrom,&
       dkpot,dpot,dvolt,dsens,drandb

    !!!$   Schalter ob transformierte Potentialwerte ausgegeben werden sollen
    logical             lkpot

    !!!$   Schalter ob Potentialwerte ausgegeben werden sollen
    logical             lpot

    !!!$   Schalter ob Spannungswerte ausgegeben werden sollen
    logical             lvolt

    !!!$   Schalter ob weiterer Datensatz modelliert werden soll
    logical             lagain

    !!!$   Schalter ob mit K-Faktor modelliert werden soll (default ohne)
    logical             wkfak

    !!!$   Schalter ob nur analytisch modelliert werden soll (default ohne)
    logical             lana

    !!!$   Indexvariablen
    integer             j,k,l,c1,info

    !!!$ counting wavenumbers
    INTEGER :: count

    character(len=256)    ftext
    INTEGER  getpid,pid,maxthreads,mythreads

    !!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    kanal = 11
    pid = getpid()
    fetxt = 'crmod.pid'
    PRINT*,'CRMod Process_ID ::',pid
    OPEN (kanal,FILE=TRIM(fetxt),STATUS='replace',err=999)
    WRITE (kanal,*)pid
    CLOSE (kanal)
    !  maxthreads = OMP_GET_MAX_THREADS()
    !  WRITE(6,"(a, i3)") " OpenMP max threads: ", maxthreads

    CALL get_git_ver(version)

    CALL tic(c1)
    !!!$   'crmod.cfg' oeffnen
    fetxt = 'crmod.cfg'
    errnr = 1
    open(12,file=TRIM(fetxt),status='old',err=999)

    !!!$   Allgemeine Parameter setzen


    !!!$   "singularity removal" ?
    !!!$   ak        lsr = .true.
    lsr = .false.

    !!!$   Art der Ruecktransformation
    !!!$   ak        swrtr = 1
    !!!$   ak        swrtr = 0

    !!!$   Sonstiges
    lkpot = .false.
    wkfak = .false.
    lana = .FALSE.
    lprior = .false. ! kann prior modell verrauschen
    ldc = .FALSE.
    lnsepri = .false.
    !!!$   ak        lkpot = .true.
    !!!$   ak        dkpot = '..\tmp\kpot.ctr'

    !!!$   Kontrollausgabe
    5 write(*,*)
    write(*,'(a20)') ' Reading Input-Files'

    !!!$   'crmod.cfg' einlesen
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

    GOTO 101

    100 BACKSPACE (12)

    101 lana = BTEST(mswitch,0)
    wkfak = BTEST(mswitch,1)
    lsr = BTEST(mswitch,2)
    lverb = BTEST(mswitch,10)

    IF (lana) PRINT*,'ANALYTICAL SOLUTION'
    IF (wkfak) PRINT*,'WITH K-FAKTOR'
    IF (lsr) THEN
        PRINT*,'WITH SINGULARITY REMOVAL'
        lana = .FALSE. ! analytical for singularity
        ! is indeed controlled by lsr and not lana..
        ! lana is only true for analytical solution only
    END IF

    IF (lverb) PRINT*,'VERBOSE OUTPUT'

    !!!$   Alles einlesen
    call relem(kanal,delem)
    if (errnr.ne.0) goto 999

    call relectr(kanal,delectr)
    if (errnr.ne.0) goto 999

    call rdatm(kanal,dstrom)
    if (errnr.ne.0) goto 999


    if (swrtr.eq.0) then
        lana = .FALSE.
        lsr = .FALSE.
        kwnanz = 1
        ALLOCATE (kwn(kwnanz),stat=errnr)
        if (errnr /= 0) then
           fetxt = 'Error memory allocation kwn'
           errnr = 94
           GOTO 999
        end if
        kwn(1) = 0d0
    else
        call rwaven()
        if (errnr.ne.0) goto 999
    end if
    print*,''
    NTHREADS = maxthreads
    !!!$ now that we know nf and kwnanz, we can adjust the OMP environment..
    if (maxthreads > 2) then ! single or double processor machines don't need scheduling..
        mythreads = MAX(kwnanz,2)
        mythreads = kwnanz
        PRINT*,'Rescheduling..'
        if ( mythreads <= maxthreads ) then ! best case,
            !!!$ the number of processors is greater or equal the assumed
            !!!$ workload
            PRINT*,'perfect match'
        else
            !!!$ is smaller than the minimum workload.. now we have to devide a bit..
            PRINT*,'less nodes than wavenumbers'
            do k = 1, INT(kwnanz/2)
                mythreads = INT(kwnanz / k) + 1
                if (mythreads < maxthreads) exit
            end do
        end if
        NTHREADS = mythreads
    end if

    nthreads = 2

    !  CALL OMP_SET_NUM_THREADS ( NTHREADS )
    ! recheck ..
    !  k = OMP_GET_MAX_THREADS()
    !  WRITE(6,'(2(a, i3),a)') " OpenMP threads: ",k,'(',maxthreads,')'
    print*
    !!!$   Startmodell belegen
    ALLOCATE (sigma(elanz),stat=errnr)
    if (errnr /= 0) then
        fetxt = 'Error memory allocation fem sigma'
        errnr = 94
        goto 999
    end if

    !!$ manz = elanz      
    manz =elanz

    call rsigma(kanal,dsigma)
    if (errnr.ne.0) goto 999

    if (lrandb2) then
        call rrandb(kanal,drandb)
        if (errnr.ne.0) goto 999
    end if

    !!!$   Ggf. Referenzleitfaehigkeit bestimmen
    call refsig()

    if (lsink) write(6,'(/A,I5,2F12.3/)')&
       'Fictious sink @ node ',nsink,sx(snr(nsink)),sy(snr(nsink))
    !!!$   Element- und Randelementbeitraege sowie Konfigurationsfaktoren zur
    !!!$   Berechnung der gemischten Randbedingung bestimmen
    call precal()
    if (errnr.ne.0) goto 999

    if (.not.lbeta) lsr=.false.

    !!!$   'a', 'hpot' und 'kpot' zuweisen
    ALLOCATE(a((mb+1)*sanz),hpot(sanz,eanz),&
       kpot(sanz,eanz,kwnanz),b(sanz),stat=errnr)
           allocate (a_mat(sanz,sanz),STAT=errnr)
           allocate (a_mat_band(3*mb+1,sanz))
           allocate(b_mat(sanz,eanz))
           allocate (ipiv(sanz),STAT=errnr)

    if (errnr.ne.0) then
        fetxt = 'allocation problem a and hpot'
        errnr = 97
        goto 999
    end if
    ALLOCATE(pot(sanz),pota(sanz),fak(sanz),stat=errnr)
    if (errnr.ne.0) then
        fetxt = 'allocation problem pot to fak'
        errnr = 97
        goto 999
    end if
    count = 0
    !  !$OMP PARALLEL DEFAULT (none) &
    !  !$OMP FIRSTPRIVATE (pota,fak,pot,a,b,fetxt) &
    !  !$OMP PRIVATE (j,l,k,a_mat,b_mat,a_mat_band,ipiv,info,a_mat_band_elec) &
    !  !$OMP SHARED (kwnanz,lverb,eanz,lsr,enr,lbeta,lrandb,lrandb2,&
    !  !$OMP  sanz,kpot,swrtr,hpot,count,lana,kg,elbg,relanz,sigma,mb)
    !  !$OMP DO
    !!!$   POTENTIALWERTE BERECHNEN
    do k=1,kwnanz
        !!!$   Kontrollausgabe
        !     !$OMP ATOMIC
        count = count + 1

        if (swrtr.eq.0) then
            write(*,'(a)')' Calculating Potentials'
        else
            write (*,'(a,t45,I4,t100,a)',ADVANCE='no')&
                   ACHAR(13)//' Calculating Potentials : Wavenumber ',count
        end if
        call pre_comp_ab(k,a_mat_band)
        do l=1,eanz
            b = cmplx(0.)
            a_mat_band_elec = a_mat_band
            b(enr(l)) = cmplx(1.)
            call comp_ab(k,a_mat_band_elec,l)
            !!!!$   Ggf. Randbedingung beruecksichtigen
            !           if (lrandb) call randb(a,b)
            !           if (lrandb2) call randb2(a,b)

            ! General Band matrix
            call zgbsv(sanz,mb,mb, 1, a_mat_band_elec, 3*mb+1, ipiv, b, sanz,info )
            if (info.ne.0) print*,'ZGBSV info:',info
            !!!$   Potentialwerte zurueckskalieren und umspeichern sowie ggf.
            !!!$   analytische Loesung addieren
            do j=1,sanz
                if (lana) then
                    kpot(j,l,k) = pota(j)
                else
                    kpot(j,l,k) = b(j)
                    if (lsr) kpot(j,l,k) = kpot(j,l,k) + pota(j)
                end if
                if (swrtr.eq.0) hpot(j,l) = kpot(j,l,k)
            end do
        end do
    end do
    !  !$OMP END PARALLEL
    print*,'done, now processing'
    !!!$   Ggf. Ruecktransformation der Potentialwerte
    if (swrtr.eq.1) call rtrafo()

    !!!$   Ggf. transformierte Potentialwerte ausgeben
    if (lkpot) then
        call wkpot(kanal,dkpot)
        if (errnr.ne.0) goto 999
    end if

    !!!$   Ggf. Potentialwerte berechnen und ausgeben
    if (lpot) then
        call bpot(kanal,dpot)
        if (errnr.ne.0) goto 999
    end if

    !!!$   Ggf. Spannungswerte berechnen und ausgeben
    !!!$   (bzw. scheinbaren Widerstandswerte)
    if (lvolt) then
        if (wkfak.and.lbeta) then
            call bkfak()
            if (errnr.ne.0) goto 999
        else
            !!!$   nur echte Spannungen ausgeben...
            kfak=1.0 
        end if

        call bvolt()
        if (errnr.ne.0) goto 999

        call wdatm(kanal,dvolt)
        if (errnr.ne.0) goto 999
    end if

    !!!$   'a' und 'hpot' freigeben
    deallocate(a,hpot,b,a_mat,b_mat,ipiv)

    !!!$   Ggf. Sensitivitaeten aller Messungen berechnen und ausgeben
    if (lsens) then
        if (manz.ne.elanz) then
            fetxt = 'manz /= elanz .. is not implemented yet'
            errnr = 50
            goto 999
        end if
        !     !$ get memory for mnr..
        ALLOCATE (mnr(elanz),stat=errnr)
        if (errnr /= 0) then
            fetxt = 'Error memory allocation mnr failed'
            errnr = 94
            goto 999
        end if
        !!!$   Modelleinteilung gemaess Elementeinteilung belegen

        do j=1,elanz
            mnr(j) = j
        end do

        !!!$   'sens' zuweisen
        ALLOCATE(sens(nanz,manz),stat=errnr)
        if (errnr.ne.0) then
            errnr = 97 
            goto 999
        end if

        call bsens()

        call wsens(kanal,dsens)
        if (errnr.ne.0) goto 999

        !!!$   'sens' freigeben
        DEALLOCATE(sens)
    end if

    !!!$   'kpot' freigeben
    DEALLOCATE(kpot,pot,pota,fak)

    !!!$   Ggf. weiteren Datensatz modellieren
    if (lagain) goto 5

    !!!$   'crmod.cfg' schliessen
    close(12)

    !!!$   Kontrollausgabe

    fetxt = 'solution time'
    CALL TOC(c1,fetxt)

    write(*,*)
    write(*,'(a)',ADVANCE='no')' Modelling completed'


    if (ALLOCATED (snr)) DEALLOCATE (snr,sx,sy)
    if (ALLOCATED (typ)) DEALLOCATE (typ,nelanz,selanz)
    if (ALLOCATED (nrel)) DEALLOCATE (nrel,rnr)

    if (ALLOCATED (strnr)) DEALLOCATE (strnr,strom,volt,sigmaa,&
       kfak,vnr)
    if (ALLOCATED (sigma)) DEALLOCATE (sigma)
    if (ALLOCATED (enr)) DEALLOCATE (enr)
    if (ALLOCATED (mnr)) DEALLOCATE (mnr)
    if (ALLOCATED (kwn)) DEALLOCATE (kwn)
    if (ALLOCATED (kwnwi)) DEALLOCATE (kwnwi)

    if (ALLOCATED (rwddc)) DEALLOCATE (rwddc) 
    if (ALLOCATED (rwndc)) DEALLOCATE (rwndc) 
    if (ALLOCATED (rwd)) DEALLOCATE (rwd) 
    if (ALLOCATED (rwn)) DEALLOCATE (rwn) 
    if (ALLOCATED (rwdnr)) DEALLOCATE (rwdnr) 

    STOP '0'

    !!!$....................................................................

    !!!$   (Fehler-) Meldung schreiben
    999 open(9,file='error.dat',status='replace')
    errflag = 2
    CALL get_error(ftext,errnr,errflag,fetxt)
    write(9,'(a80,i3,i1)') fetxt,errnr,errflag
    write(9,*)ftext
    close(9)
    stop '-1'

    1001 open(9,file='error.dat',status='replace')
    errnr   = 2
    errflag = 2
    CALL get_error(ftext,errnr,errflag,fetxt)
    write(9,'(a80,i3,i1)') fetxt,errnr,errflag
    write(9,*)ftext
    close(9)
    stop '-2'

end program fem

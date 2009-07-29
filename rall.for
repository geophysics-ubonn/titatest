        subroutine rall(kanal,delem,delectr,dstrom,drandb,
cdiff-     1                  dsigma,dvolt,dsens,dstart,lsens,lagain)
cdiff+<
     1                  dsigma,dvolt,dsens,dstart,dd0,dm0,dfm0,
     1                  lsens,lagain)
cdiff+>
    
c Unterprogramm zum Einlesen der benoetigten Variablen.

c Andreas Kemna                                            01-Mar-1995
c                                       Letzte Aenderung   20-Aug-2007

c.....................................................................

        INCLUDE 'parmax.fin'
        INCLUDE 'err.fin'
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
        INCLUDE 'randb.fin'

c.....................................................................

c EIN-/AUSGABEPARAMETER:

c Kanalnummer
        integer         * 4     kanal

c Dateinamen
        character       * 80    delem,
     1                          delectr,
     1                          dstrom,
     1                          dsigma,
     1                          dvolt,
     1                          dsens,
     1                          dstart,
cdiff+<
     1                          dd0,
     1                          dm0,
     1                          dfm0,
cdiff+>
     1                          drandb

c Schalter ob Summe der Sensitivitaeten aller Messungen ausgegeben
c werden soll
        logical         * 4     lsens

c Schalter ob weiterer Datensatz invertiert werden soll
        logical         * 4     lagain

c.....................................................................

c PROGRAMMINTERNE PARAMETER:

c Indexvariable
        integer         * 4     i

c Pi
        real            * 8     pi

cdiff+<
        real            * 8     dum(nmax),dum2(nmax),dum3
        integer         * 4     ic(nmax),ip(nmax),idum(nmax),nanz0,j,j0
cdiff+>

cak Inga
        integer         * 4     elec1,elec2,
     1                          elec3,elec4

c.....................................................................

        pi = dacos(-1d0)

c 'crtomo.cfg' EINLESEN
        fetxt = 'crtomo.cfg'
        errnr = 3
        read(12,*,end=1001,err=999)
        read(12,'(a80)',end=1001,err=999) delem
        read(12,'(a80)',end=1001,err=999) delectr
        read(12,'(a80)',end=1001,err=999) dstrom
        read(12,'(a60)',end=1001,err=999) ramd
cdiff+<
        read(12,*,end=1001,err=999) ldiff
        read(12,'(a80)',end=1001,err=999) dd0
        read(12,'(a80)',end=1001,err=999) dm0
        read(12,'(a80)',end=1001,err=999) dfm0
cdiff+>
        read(12,*,end=1001,err=999)
        read(12,*,end=1001,err=999) nx
        read(12,*,end=1001,err=999) nz
        read(12,*,end=1001,err=999) alfx
        read(12,*,end=1001,err=999) alfz
        read(12,*,end=1001,err=999) itmax
cak        read(12,*,end=1001,err=999) nrmsdm
        read(12,*,end=1001,err=999) ldc
cak        read(12,*,end=1001,err=999) lsr
        read(12,*,end=1001,err=999) lrobust
cak        read(12,*,end=1001,err=999) lpol
        read(12,*,end=1001,err=999) lfphai
cak        read(12,*,end=1001,err=999) lindiv
        read(12,*,end=1001,err=999) stabw0
        read(12,*,end=1001,err=999) stabm0
        read(12,*,end=1001,err=999) stabpA1
        read(12,*,end=1001,err=999) stabpB
        read(12,*,end=1001,err=999) stabpA2
        read(12,*,end=1001,err=999) stabp0
        read(12,*,end=1001,err=999) lrho0
        read(12,*,end=1001,err=999) bet0
        read(12,*,end=1001,err=999) pha0
        read(12,*,end=1001,err=999) lagain
        read(12,*,end=1001,err=999) swrtr
        read(12,*,end=1001,err=999) lsink
        read(12,*,end=1001,err=999) nsink
        read(12,*,end=1001,err=999) lrandb2
        read(12,'(a80)',end=1001,err=999) drandb

cro        lsr    = .false.
cro        lpol   = .true.
cro        lfphai = .true.
cro        lrho0  = .false.

cakERT2003
cak        ldc    = .true.
cak        lsr    = .false.
cak        lpol   = .false.

        nrmsdm = 1d0
        lsr    = .false.
        lpol   = .false.
        lindiv = .false.
        ltri   = 1

c Ggf. Fehlermeldungen
        if (nx.lt.2.or.nz.lt.2) then
            fetxt = ' '
            errnr = 89
            goto 999
        else if (alfx.le.0d0.or.alfz.le.0d0) then
            fetxt = ' '
            errnr = 96
            goto 999
        else if (itmax.lt.1.or.itmax.ge.100) then
            fetxt = ' '
            errnr = 61
            goto 999
        else if (nrmsdm.lt.1d-12) then
            fetxt = ' '
            errnr = 62
            goto 999
c        else if (nlam.lt.0) then
c            fetxt = ' '
c            errnr = 83
c            goto 999
c        else if (fstart.gt.1d0.or.fstop.gt.1d0.or.
c     1           fstart.le.0d0.or.fstop.le.0d0.or.
c     1           fstart.gt.fstop) then
c            fetxt = ' '
c            errnr = 98
c            goto 999
        else if (stabw0.le.0d0.or.stabm0.lt.0d0) then
            fetxt = ' '
            errnr = 104
            goto 999
        else if (.not.ldc.and.lfphai.and.
     1           (stabp0.le.0d0.or.stabpA2.lt.0d0)) then
            fetxt = ' '
            errnr = 105
            goto 999
        else if (lrho0.and.(bet0.le.0d0.or.
     1                      (.not.ldc.and.dabs(pha0).gt.1d3*pi))) then
            fetxt = ' '
            errnr = 91
            goto 999
c        else if (mqrms.lt.0d0.or.mqrms.ge.1d0) then
c            fetxt = ' '
c            errnr = 64
c            goto 999
c        else if (lrobust.and.l1min.lt.1d0) then
c            fetxt = ' '
c            errnr = 90
c            goto 999
        end if
        
c FIXED PARAMETER
c Slash
        slash = '/'

cak
c 'dstart'
        lstart = .false.
        dstart = ' '
cak        lstart = .true.
cak        dstart = '..\..\strasbrg\9610\plane45\mod\rho0.dat'

c "Force negative phase" ?
csandra        lphi0 = .true.
        lphi0 = .true.
cak        lphi0 = .false.

c "ratio-dataset" ?
        lratio = .false.
cak        lratio = .true.

c Minimale "L1-ratio" (Grenze der "robust inversion")
        l1min = 1d0
cak        l1min = 1.2d0

c Ueberdeckung schreiben ?
c2d        lsens = .true.
        lsens = .false.

c Art der Ruecktransformation
cak        swrtr = 1

c Minimaler Quotient zweier aufeinanderfolgender Daten-RMS-Werte
cak Default
cak        mqrms = 1d-2
cak ERT 2002-2003 synth
        mqrms = 2d-2

cak        mqrms = 2d-2
cak Tank
cak        mqrms = 2d-2
cak MMAJ
cak        mqrms = 5d-2

c CG-Epsilon
        eps = 1d-4

c Mindest-step-length
        stpmin = 1d-3

c Regularisierungsparameter
cak Default
        nlam   = 20

cak Default
        fstart = 0.5d0
        fstop  = 0.9d0

cak MMAJ
cak        fstart = 0.2d0
cak        fstop  = 0.8d0
cak Strasbrg/Werne/Grimberg
cak        fstart = 0.5d0
cak        fstop  = 0.8d0

c Sonstiges
        if (lratio) then
            lrho0  = .true.
            lstart = .false.
            lphi0  = .false.
            lpol   = .false.
        end if

cdiff-        if (lstart) lrho0=.false.
cdiff+<
        if (lstart.or.ldiff) lrho0=.false.
cak
        if (ldiff) then
            ldc  = .true.
            lpol = .false.
        end if
cdiff+>
cak        if (ldc.or.stabp0.ge.stabw0) lfphai=.false.
        if (ldc) lfphai=.false.

c Dateien
        lnramd = index(ramd,' ')-1
        dsigma = ramd(1:lnramd)//slash(1:1)//'rho.dat'
        dvolt  = ramd(1:lnramd)//slash(1:1)//'volt.dat'
        dsens  = ramd(1:lnramd)//slash(1:1)//'sens.dat'

c Elementeinteilung einlesen
        call relem(kanal,delem)
        if (errnr.ne.0) goto 999
        
        IF (ltri==1) CALL bnachbar
c Modelleinteilung gemaess Elementeinteilung belegen
        manz = nx*nz

        if (manz.ne.elanz) then
            fetxt = ' '
            errnr = 50
            goto 999
        else if (manz.gt.mmax) then
            fetxt = ' '
            errnr = 63
            goto 999
        end if

        do i=1,elanz
            mnr(i) = i
        end do

c Maximale Anzahl an CG-steps setzen
cak        ncgmax = manz
        ncgmax = manz/10

c Elektrodenverteilung und Daten einlesen sowie Wellenzahlwerte
c bestimmen
        call relectr(kanal,delectr)
        if (errnr.ne.0) goto 999

        call rdati(kanal,dstrom)
        if (errnr.ne.0) goto 999

        if (swrtr.eq.0) then
            lsr    = .false.
            kwnanz = 1
            kwn(1) = 0d0
        else
            call rwaven()
            if (errnr.ne.0) goto 999
        end if

c read boundary values
        if (lrandb2) then
            call rrandb(kanal,drandb)
            if (errnr.ne.0) goto 999
        end if

cdiff+<
        if (ldiff) then
            open(kanal,file=dd0,status='old')
            read(kanal,*) nanz0
            do j=1,nanz0

                read(kanal,*) ic(j),ip(j),dum(j)
cak Inga
cak                read(kanal,*) elec1,elec2,elec3,elec4,dum(j)
cak                ic(j) = elec1*10000 + elec2
cak                ip(j) = elec3*10000 + elec4

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
10          i  = i+1
                j = j0
20              j = j+1
                if (strnr(i).eq.ic(j).and.vnr(i).eq.ip(j).and.
     1              idum(j).eq.1) then
c nur falls jede Messkonfiguration nur einmal!
c                    j0     = j
                    d0(i)  = dcmplx(-dlog(dum(j)),0d0)
                    fm0(i) = dcmplx(-dlog(dum2(j)),0d0)
                else if (j.lt.nanz0) then
                    goto 20
                else
                    errnr = 1
                    fetxt = ramd(1:lnramd)//slash(1:1)//'run.ctr'
                    open(10,file=fetxt,status='unknown',err=999)
1                   read(10,*,end=2)
                        goto 1
2                   backspace(10)
                    errnr = 4
                    write(10,'(i7,1x,i7,a12)',err=999)
     1                       strnr(i),vnr(i),' : discarded'
                    close(10)

                    nanz = nanz-1
                    do j=i,nanz
                        strnr(j) = strnr(j+1)
                        vnr(j)   = vnr(j+1)
                        dat(j)   = dat(j+1)
                        wmatd(j) = wmatd(j+1)
                        if (lfphai) wmatdp(j)=wmatdp(j+1)
c nicht notwendig, da Werte eh alle 1
c                        wdfak(j) = wdfak(j+1)
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
        end if
cdiff+>

        errnr = 0
        return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c Fehlermeldungen

999     return

1001    errnr = 2
        return

        end

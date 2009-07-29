        subroutine rrandb(kanal,datei)

c Unterprogramm zum Einlesen der Randwerte aus 'datei'.

c Andreas Kemna                                            12-Feb-1993
c                                       Letzte Aenderung   15-Jul-2007

c.....................................................................

        INCLUDE 'parmax.fin'
        INCLUDE 'err.fin'
        INCLUDE 'elem.fin'
        INCLUDE 'fem.fin'
        INCLUDE 'randb.fin'

c.....................................................................

c EIN-/AUSGABEPARAMETER:

c Kanalnummer
        integer         * 4     kanal

c Datei
        character       * 80    datei

c.....................................................................

c PROGRAMMINTERNE PARAMETER:

c Indexvariable
        integer         * 4     i

c.....................................................................

c 'datei' oeffnen
        fetxt = datei

        errnr = 1
        open(kanal,file=fetxt,status='old',err=999)
        errnr = 3

c DC CASE
        if (ldc) then

c Anzahl der Randwerte einlesen (Dirichlet)
        read(kanal,*,end=1001,err=1000) rwdanz

c Ggf. Fehlermeldung
        if (rwdanz.gt.rwdmax) then
            fetxt = ' '
            errnr = 13
            goto 1000
        end if

c Knotennummern der Randwerte sowie Randwerte einlesen (Dirichlet)
        do i=1,rwdanz
            read(kanal,*,end=1001,err=1000) rwdnr(i),rwddc(i)

c Ggf. Fehlermeldung
            if (rwdnr(i).gt.sanz) then
                fetxt = ' '
                errnr = 30
                goto 1000
            end if
        end do

c Anzahl der Randwerte einlesen (Neumann)
        read(kanal,*,end=1001,err=1000) rwnanz

c Ggf. Fehlermeldung
        if (rwnanz.gt.relmax) then
            fetxt = ' '
            errnr = 33
            goto 1000
        end if

c Randwerte einlesen (Neumann)
        if (rwnanz.gt.0)
     1      read(kanal,*,end=1001,err=1000) (rwndc(i),i=1,rwnanz)
        else

c COMPLEX CASE

c Anzahl der Randwerte einlesen (Dirichlet)
        read(kanal,*,end=1001,err=1000) rwdanz

c Ggf. Fehlermeldung
        if (rwdanz.gt.rwdmax) then
            fetxt = ' '
            errnr = 13
            goto 1000
        end if

c Knotennummern der Randwerte sowie Randwerte einlesen (Dirichlet)
        do i=1,rwdanz
            read(kanal,*,end=1001,err=1000) rwdnr(i),rwd(i)

c Ggf. Fehlermeldung
            if (rwdnr(i).gt.sanz) then
                fetxt = ' '
                errnr = 30
                goto 1000
            end if
        end do

c Anzahl der Randwerte einlesen (Neumann)
        read(kanal,*,end=1001,err=1000) rwnanz

c Ggf. Fehlermeldung
        if (rwnanz.gt.relmax) then
            fetxt = ' '
            errnr = 33
            goto 1000
        end if

c Randwerte einlesen (Neumann)
        if (rwnanz.gt.0)
     1      read(kanal,*,end=1001,err=1000) (rwn(i),i=1,rwnanz)

c 'datei' schliessen
        end if

        close(kanal)

        errnr = 0
        return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c Fehlermeldungen

999     return

1000    close(kanal)
        return

1001    close(kanal)
        errnr = 2
        return

        end

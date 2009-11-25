        subroutine relectr(kanal,datei)
 
c Unterprogramm zum Einlesen der Elektrodenverteilung aus 'datei'.

c Andreas Kemna                                            11-Oct-1993
c                                       Letzte Aenderung   24-Oct-1996

c.....................................................................

        INCLUDE 'parmax.fin'
        INCLUDE 'err.fin'
        INCLUDE 'elem.fin'
        INCLUDE 'electr.fin'

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

c Anzahl der Elektroden einlesen
        read(kanal,*,end=1001,err=1000) eanz
c Ggf. Fehlermeldung
        if (eanz.gt.emax) then
            fetxt = ' '
            errnr = 12
            goto 1000
        end if

c Knotennummern der Elektroden einlesen
        do i=1,eanz
            read(kanal,*,end=1001,err=1000) enr(i)

c Ggf. Fehlermeldung
            if (enr(i).gt.sanz) then
                fetxt = ' '
                errnr = 29
                goto 1000
            end if
        end do

c 'datei' schliessen
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

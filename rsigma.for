        subroutine rsigma(kanal,datei)

c Unterprogramm zum Einlesen der Widerstandsverteilung aus 'datei'.

c Andreas Kemna                                            20-Dec-1993
c                                       Letzte Aenderung   07-Nov-1997
     
c.....................................................................

        INCLUDE 'parmax.fin'
        INCLUDE 'err.fin'
        INCLUDE 'elem.fin'
        INCLUDE 'sigma.fin'

c.....................................................................

c EIN-/AUSGABEPARAMETER:

c Kanalnummer
        integer         * 4     kanal

c Datei
        character       * 80    datei

c.....................................................................

c PROGRAMMINTERNE PARAMETER:

c Hilfsvariablen
        integer         * 4     i,idum
        real            * 8     bet,pha

c.....................................................................

c 'datei' oeffnen
        fetxt = datei

        errnr = 1
        open(kanal,file=fetxt,status='old',err=999)

        errnr = 3

c Anzahl der Elemente (ohne Randelemente) einlesen
        read(kanal,*,end=1001,err=1000) idum

c Ggf. Fehlermeldung
        if (idum.ne.elanz) then
            fetxt = ' '
            errnr = 47
            goto 1000
        end if

c Betrag und Phase (in mrad) des komplexen Widerstandes einlesen
        do i=1,elanz
            read(kanal,*,end=1001,err=1000) bet,pha

c Ggf. Fehlermeldung
            if (bet.lt.1d-12) then
                fetxt = ' '
                errnr = 11
                goto 999
            end if

c Komplexe Leitfaehigkeit bestimmen
            pha      = 1d-3*pha
            sigma(i) = dcmplx(dcos(pha)/bet,-dsin(pha)/bet)
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

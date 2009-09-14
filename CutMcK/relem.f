        subroutine relem(kanal,datei)

c Unterprogramm zum Einlesen der FEM-Parameter aus 'datei'.

c Andreas Kemna                                            11-Oct-1993
c                                       Letzte Aenderung   24-Oct-1996

c.....................................................................

        INCLUDE 'parmax.fin'
        INCLUDE 'err.fin'
        INCLUDE 'elem.fin'
                 
c.....................................................................

c EIN-/AUSGABEPARAMETER:

c Kanalnummer
        integer         * 4     kanal

c Datei
        character       * 80    datei

c.....................................................................

c PROGRAMMINTERNE PARAMETER:

c Indexvariablen
        integer         * 4     i,j,k

c Hilfsvariable
        integer         * 4     idum

c.....................................................................

c 'datei' oeffnen
        fetxt = datei

        errnr = 1
        open(kanal,file=fetxt,status='old',err=999)

        errnr = 3

c Anzahl der Knoten (bzw. Knotenvariablen), Anzahl der Elementtypen
c sowie Bandbreite der Gesamtsteifigkeitsmatrix einlesen
        read(kanal,*,end=1001,err=1000) sanz,typanz,mb

c Ggf. Fehlermeldungen
        if (sanz.gt.smax) then
            fetxt = ' '
            errnr = 5
            goto 1000
        else if (typanz.gt.typmax) then
            fetxt = ' '
            errnr = 6
            goto 1000
        else if (mb.gt.mbmax) then
            fetxt = ' '
            errnr = 7
            goto 1000
        end if

c Elementtypen, Anzahl der Elemente eines bestimmten Typs sowie
c Anzahl der Knoten in einem Elementtyp einlesen
        read(kanal,*,end=1001,err=1000)
     1             (typ(i),nelanz(i),selanz(i),i=1,typanz)

c Anzahl der Elemente (ohne Randelemente) und Anzahl der Randelemente
c bestimmen
        relanz = 0
        elanz  = 0
        do i=1,typanz
            if (typ(i).gt.10) then
                relanz = relanz + nelanz(i)
            else
                elanz  = elanz  + nelanz(i)
            end if

c Ggf. Fehlermeldung
            if (selanz(i).gt.selmax) then
                fetxt = ' '
                errnr = 8
                goto 1000
            end if
        end do
c Ggf. Fehlermeldungen
        if (elanz.gt.elmax) then
            fetxt = ' '
            errnr = 9
            goto 1000
        else if (relanz.gt.relmax) then
            fetxt = ' '
            errnr = 10
            goto 1000
        end if

        PRINT*,'check dis out ',4+sanz+elanz+relanz*2

c Zeiger auf Koordinaten, x-Koordinaten sowie y-Koordinaten der Knoten
c einlesen
        print*,'knoten::',sanz
        read(kanal,*,end=1001,err=1000) (snr(i),sx(i),sy(i),i=1,sanz)
        i = 0;j = 0
        DO k=1,sanz
           IF (sx(snr(k))==sx(snr(1))) j = j+1 !counts vertical grid nodes
           IF (sx(snr(k))==sx(snr(1))) i = i+1 !counts horizontal grid nodes
        END DO
        WRITE (*,'(3(a,I6))')'Counted nx ',i,' ny',j,' nxy=',i*j
c Knotennummern der Elemente einlesen
        print*,'post knoten position',4+sanz
        idum = 0
        print*,'elemente::',elanz,relanz
        do i=1,typanz
            do j=1,nelanz(i)
                read(kanal,*,end=1001,err=1000)
     1                     (nrel(idum+j,k),k=1,selanz(i))
            end do
            idum = idum + nelanz(i)
        end do
        print*,'post elemente position',4+sanz+elanz+relanz
        print*,'randpunkte::',relanz
c Zeiger auf Werte der Randelemente einlesen
        read(kanal,*,end=1001,err=1000) (rnr(i),i=1,relanz)
        print*,'post randpunkte position',4+sanz+elanz+relanz*2
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

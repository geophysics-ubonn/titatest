        subroutine bpot(kanal,datei)

c Unterprogramm zur Berechnung der Potentialwerte
c (beachte: Potentialwerte wurden fuer Einheitsstrom berechnet).

c Andreas Kemna                                            17-Aug-1994
c                                       Letzte Aenderung   24-Jun-1997

c.....................................................................

        USE alloci

        INCLUDE 'parmax.fin'
        INCLUDE 'err.fin'
        INCLUDE 'elem.fin'
        INCLUDE 'dat.fin'
        INCLUDE 'fem.fin'

c.....................................................................

c EIN-/AUSGABEPARAMETER:

c Kanalnummer
        integer         * 4     kanal

c Datei
        character       * 80    datei

c.....................................................................

c PROGRAMMINTERNE PARAMETER:

c Indexvariablen
        integer         * 4     i,j

c Elektrodennummern
        integer         * 4     elec1,elec2

c.....................................................................

        do i=1,nanz

c Stromelektroden bestimmen
            elec1 = mod(strnr(i),10000)
            elec2 = (strnr(i)-elec1)/10000

            do j=1,sanz

c (beachte: Faktoren '.../2d0' (-> Potentialwerte fuer Einheitsstrom)
c  und '...*2d0' (-> Ruecktransformation) kuerzen sich weg !)
                if (elec1.eq.0) then
                    pot(j) = hpot(j,elec2) * dcmplx(strom(i))
                else if (elec2.eq.0) then
                    pot(j) = -hpot(j,elec1) * dcmplx(strom(i))
                else
                    pot(j) = (hpot(j,elec2)-hpot(j,elec1))
     1                       * dcmplx(strom(i))
                end if
            end do

c Potentialwerte ausgeben
            call wpot(kanal,datei,i)
            if (errnr.ne.0) goto 1000
        end do

        errnr = 0
        return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c Fehlermeldungen

1000    return

        end

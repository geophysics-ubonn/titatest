        real * 4 function chareal(txt,ltxt)

c Die Funktion wandelt einen String in einen Real-Wert;
c bei Real-Werten muss die Laenge incl. '.' angegeben werden.

c Andreas Kemna                                            03-Feb-1993
c                                       Letzte Aenderung   24-Oct-1996
           
c.....................................................................

        character       * (*)   txt
        integer         * 4     ltxt
        integer         * 4     i,ih1,pos,j
        real            * 4     h2
        logical         * 4     neg

c.....................................................................

        h2  = 0.
        pos = index(txt(1:ltxt),'.')-1

c Vorzeichen bestimmen
        if (index(txt(1:ltxt),'-').gt.0) then
            neg = .true.
        else
            neg = .false.
        end if

c Umrechnung von Integer-Werten ermoeglichen
        if (pos.eq.-1) pos=ltxt

        j = 0

c Vorkommastellen
        do 10 i=pos,1,-1
            ih1 = ichar(txt(i:i))

            if (ih1.ge.48.and.ih1.le.57) then
                h2 = real(ih1-48)*10.**j + h2
                j  = j+1
            end if
10      continue

c Dezimalstellen
        if (pos.ne.ltxt) then
            j = -1

            do 20 i=pos+2,ltxt
                ih1 = ichar(txt(i:i))

                if (ih1.ge.48.and.ih1.le.57) then
                    h2 = real(ih1-48)*10.**j + h2
                    j  = j-1
                end if

20          continue
        end if

        if (neg) then
            chareal = -h2
        else
            chareal = h2
        end if

        return
        end
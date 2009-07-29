        character * 12 function intcha(idum,lidum)

c Die Function wandelt einen Integer 'idum' in einen String der
c Laenge 'lidum'.

c Andreas Kemna                                            22-Jan-1993
c                                       Letzte Aenderung   24-Oct-1996

c.....................................................................

        integer         * 4     idum
        integer         * 4     lidum

        integer         * 4     i,ih1
        character       * 5     form

c.....................................................................

c Integer in String schreiben
        if (lidum.lt.10) then
            write(form,'(a2,i1,a1)') '(i',lidum,')'
        else
            write(form,'(a2,i2,a1)') '(i',lidum,')'
        end if

        write(intcha,form) idum

c Moegliche Blanks mit Nullen ueberschreiben
        do 10 i=1,lidum
            ih1 = ichar(intcha(i:i))

            if (ih1.eq.32) then
                intcha(i:i) = '0'
            end if
10      continue

        return
        end
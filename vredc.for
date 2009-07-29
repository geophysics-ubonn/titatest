        subroutine vredc()

c Fuehrt das Vorwaerts- und Rueckwaertseinsetzen mit der Cholesky-Links-
c dreiecksmatrix aus;
c 'bdc' bleibt unveraendert, 'pot' ist Loesungsvektor.

c ( Vgl. Subroutine 'VRBNDN' in Schwarz (1991) )
  
c Andreas Kemna                                            11-Oct-1993
c                                       Letzte Aenderung   14-Nov-1997

c.....................................................................

        USE alloci

        INCLUDE 'parmax.fin'
        INCLUDE 'elem.fin'
        INCLUDE 'fem.fin'

c.....................................................................

c PROGRAMMINTERNE PARAMETER:

c Hilfsvariablen
        real            * 8     potdc(smax)
        integer         * 4     idi,i0
        integer         * 4     m1,jlow
        real            * 8     s

c Indexvariablen
        integer         * 4     i,j

c.....................................................................

        m1 = mb+1

        do 20 i=1,sanz
            idi  = i*m1
            s    = bdc(i)
            i0   = idi-i
            jlow = max0(1,i-mb)

            do 10 j=jlow,i-1
                s = s - adc(i0+j)*potdc(j)
10          continue

            potdc(i) = s / adc(idi)
20      continue

        do 40 i=sanz,1,-1
            potdc(i) = - potdc(i) / adc(idi)

            jlow = max0(1,i-mb)
            i0   = idi-i

            do 30 j=jlow,i-1
                potdc(j) = potdc(j) + adc(i0+j)*potdc(i)
30          continue

            idi = idi-m1
40      continue

        do i=1,sanz
            pot(i) = dcmplx(potdc(i))
        end do

        return
        end

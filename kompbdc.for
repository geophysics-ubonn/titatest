        subroutine kompbdc(nelec)

c Unterprogramm zur Kompilation des Konstanten- bzw. Stromvektors 'bdc'
c fuer Einheitsstrom.

c Andreas Kemna                                            11-Oct-1993
c                                       Letzte Aenderung   16-Jul-2007

c.....................................................................

        INCLUDE 'parmax.fin'
        INCLUDE 'elem.fin'
        INCLUDE 'electr.fin'
        INCLUDE 'fem.fin'
           
c.....................................................................

c EIN-/AUSGABEPARAMETER:

c Aktuelle Elektrodennummer
        integer         * 4     nelec

c.....................................................................

c PROGRAMMINTERNE PARAMETER:

c Indexvariable
        integer         * 4     i

c.....................................................................

c Konstantenvektor auf Null setzen
        do i=1,sanz
            bdc(i) = 0d0
        end do

c Aufbau des Konstanten- bzw. Stromvektors mit Skalierung
c ( A * x + b = 0 )
        bdc(enr(nelec)) = -fak(enr(nelec))

cakc BAW-Tank
cak        bdc(211) = fak(211)
cakc Model EGS2003
cak        bdc(1683) = fak(1683)
cakc Lysimeter hor_elem\normal
cak        bdc(129) = fak(129)
cakc Lysimeter hor_elem\fine
cak        bdc(497) = fak(497)
cakc Simple Tucson Model
cak        bdc(431) = fak(431)
cakc TU Berlin Mesokosmos
cak        bdc(201) = fak(201)
cakc Andy
cak        bdc(2508) = fak(2508)
cakc Sandra (ele?_anom)
cak        bdc(497) = fak(497)

        if (lsink) bdc(nsink) = fak(nsink)

        return
        end

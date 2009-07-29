        subroutine kompb(nelec)

c Unterprogramm zur Kompilation des Konstanten- bzw. Stromvektors 'b'
c fuer Einheitsstrom.

c Andreas Kemna                                            11-Oct-1993
c                                       Letzte Aenderung   15-Jul-2007

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
            b(i) = dcmplx(0d0)
        end do

c Aufbau des Konstanten- bzw. Stromvektors mit Skalierung
c ( A * x + b = 0 )
        b(enr(nelec)) = dcmplx(-fak(enr(nelec)))

cakc BAW-Tank
cak        b(211) = dcmplx(fak(211))
cakc Model EGS2003
cak        b(1683) = dcmplx(fak(1683))
cakc Lysimeter hor_elem\normal
cak        b(129) = dcmplx(fak(129))
cakc Lysimeter hor_elem\fine
cak        b(497) = dcmplx(fak(497))
cakc Simple Tucson Model
cak        b(431) = dcmplx(fak(431))
cakc TU Berlin Mesokosmos
cak        b(201) = dcmplx(fak(201))
cakc Andy
cak        b(2508) = dcmplx(fak(2508))
cakc Sandra (ele?_anom)
cak        b(497) = dcmplx(fak(497))
cakc Adrian (Tank)
cak        b(1660) = dcmplx(fak(1660))

        if (lsink) b(nsink) = dcmplx(fak(nsink))

        return
        end

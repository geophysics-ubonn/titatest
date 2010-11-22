      subroutine kompb(nelec)

!!!$     Unterprogramm zur Kompilation des Konstanten- bzw. Stromvektors 'b'
!!!$     fuer Einheitsstrom.

!!!$     Andreas Kemna                                            11-Oct-1993
!!!$     Letzte Aenderung   15-Jul-2007

!!!$.....................................................................

      USE femmod
      USE electrmod
      USE elemmod

      IMPLICIT none

!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Aktuelle Elektrodennummer
      integer         * 4     nelec

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Indexvariable
      integer         * 4     i

!!!$.....................................................................

!!!$     Konstantenvektor auf Null setzen
      do i=1,sanz
         b(i) = dcmplx(0d0)
      end do

!!!$     Aufbau des Konstanten- bzw. Stromvektors mit Skalierung
!!!$     ( A * x + b = 0 )
      b(enr(nelec)) = dcmplx(-fak(enr(nelec)))

!!!$     akc BAW-Tank
!!!$     ak        b(211) = dcmplx(fak(211))
!!!$     akc Model EGS2003
!!!$     ak        b(1683) = dcmplx(fak(1683))
!!!$     akc Lysimeter hor_elem\normal
!!!$     ak        b(129) = dcmplx(fak(129))
!!!$     akc Lysimeter hor_elem\fine
!!!$     ak        b(497) = dcmplx(fak(497))
!!!$     akc Simple Tucson Model
!!!$     ak        b(431) = dcmplx(fak(431))
!!!$     akc TU Berlin Mesokosmos
!!!$     ak        b(201) = dcmplx(fak(201))
!!!$     akc Andy
!!!$     ak        b(2508) = dcmplx(fak(2508))
!!!$     akc Sandra (ele?_anom)
!!!$     ak        b(497) = dcmplx(fak(497))
!!!$     akc Adrian (Tank)
!!!$     ak        b(1660) = dcmplx(fak(1660))

      if (lsink) b(nsink) = dcmplx(fak(nsink))

      return
      end

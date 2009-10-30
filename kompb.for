      subroutine kompb(nelec)

c     Unterprogramm zur Kompilation des Konstanten- bzw. Stromvektors 'b'
c     fuer Einheitsstrom.

c     Andreas Kemna                                            11-Oct-1993
c     Letzte Aenderung   15-Jul-2007

c.....................................................................

      IMPLICIT none
      INCLUDE 'parmax.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'electr.fin'
      INCLUDE 'fem.fin'
      
c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Aktuelle Elektrodennummer
      integer         * 4     nelec

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Indexvariable
      integer         * 4     i

c.....................................................................

c     Konstantenvektor auf Null setzen
      do i=1,sanz
         b(i) = dcmplx(0d0)
      end do

c     Aufbau des Konstanten- bzw. Stromvektors mit Skalierung
c     ( A * x + b = 0 )
      b(enr(nelec)) = dcmplx(-fak(enr(nelec)))

c     akc BAW-Tank
c     ak        b(211) = dcmplx(fak(211))
c     akc Model EGS2003
c     ak        b(1683) = dcmplx(fak(1683))
c     akc Lysimeter hor_elem\normal
c     ak        b(129) = dcmplx(fak(129))
c     akc Lysimeter hor_elem\fine
c     ak        b(497) = dcmplx(fak(497))
c     akc Simple Tucson Model
c     ak        b(431) = dcmplx(fak(431))
c     akc TU Berlin Mesokosmos
c     ak        b(201) = dcmplx(fak(201))
c     akc Andy
c     ak        b(2508) = dcmplx(fak(2508))
c     akc Sandra (ele?_anom)
c     ak        b(497) = dcmplx(fak(497))
c     akc Adrian (Tank)
c     ak        b(1660) = dcmplx(fak(1660))

      if (lsink) b(nsink) = dcmplx(fak(nsink))

      return
      end

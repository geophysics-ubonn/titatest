      subroutine kompbdc(nelec)

c     Unterprogramm zur Kompilation des Konstanten- bzw. Stromvektors 'bdc'
c     fuer Einheitsstrom.

c     Andreas Kemna                                            11-Oct-1993
c     Letzte Aenderung   16-Jul-2007

c.....................................................................

      USE femmod
      USE electrmod
      USE elemmod

      IMPLICIT none

      INCLUDE 'parmax.fin'
      
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
         bdc(i) = 0d0
      end do

c     Aufbau des Konstanten- bzw. Stromvektors mit Skalierung
c     ( A * x + b = 0 )
      bdc(enr(nelec)) = -fak(enr(nelec))

c     akc BAW-Tank
c     ak        bdc(211) = fak(211)
c     akc Model EGS2003
c     ak        bdc(1683) = fak(1683)
c     akc Lysimeter hor_elem\normal
c     ak        bdc(129) = fak(129)
c     akc Lysimeter hor_elem\fine
c     ak        bdc(497) = fak(497)
c     akc Simple Tucson Model
c     ak        bdc(431) = fak(431)
c     akc TU Berlin Mesokosmos
c     ak        bdc(201) = fak(201)
c     akc Andy
c     ak        bdc(2508) = fak(2508)
c     akc Sandra (ele?_anom)
c     ak        bdc(497) = fak(497)

      if (lsink) bdc(nsink) = fak(nsink)

      return
      end

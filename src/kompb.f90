subroutine kompb(nelec,b_komp,fak_komp)

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
!!$ Berechnete Potentialwerte (bzw. Loesungsverktor)
  COMPLEX (prec),DIMENSION(*):: b_komp
!!$ Skalirerungsfaktor
  REAL (prec),DIMENSION(*)  :: fak_komp

!!!$     Aktuelle Elektrodennummer
  INTEGER (KIND = 4) ::     nelec

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Indexvariable
  INTEGER (KIND = 4) ::     i

!!!$.....................................................................

!!!$     Konstantenvektor auf Null setzen
!!$  b_komp = 0.
  do i=1,sanz
     b_komp(i) = CMPLX(0d0)
  end do

!!!$     Aufbau des Konstanten- bzw. Stromvektors mit Skalierung
!!!$     ( A * x + b = 0 )
  b_komp(enr(nelec)) = CMPLX(-fak_komp(enr(nelec)))

!!!$     akc BAW-Tank
!!!$     ak        b_komp(211) = CMPLX(fak(211))
!!!$     akc Model EGS2003
!!!$     ak        b_komp(1683) = CMPLX(fak(1683))
!!!$     akc Lysimeter hor_elem\normal
!!!$     ak        b_komp(129) = CMPLX(fak(129))
!!!$     akc Lysimeter hor_elem\fine
!!!$     ak        b_komp(497) = CMPLX(fak(497))
!!!$     akc Simple Tucson Model
!!!$     ak        b_komp(431) = CMPLX(fak(431))
!!!$     akc TU Berlin Mesokosmos
!!!$     ak        b_komp(201) = CMPLX(fak(201))
!!!$     akc Andy
!!!$     ak        b_komp(2508) = CMPLX(fak(2508))
!!!$     akc Sandra (ele?_anom)
!!!$     ak        b_komp(497) = CMPLX(fak(497))
!!!$     akc Adrian (Tank)
!!!$     ak        b_komp(1660) = CMPLX(fak(1660))

  if (lsink) b_komp(nsink) = CMPLX(fak_komp(nsink))

  return
end subroutine kompb

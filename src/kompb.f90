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
!!$ Berechnete Potentialwerte (bzw. Loesungsverktor)
!  COMPLEX (prec),DIMENSION(*,*):: b_komp
!!$ Skalirerungsfaktor
!  REAL (prec),DIMENSION(*)  :: fak_komp

!!!$     Aktuelle Elektrodennummer
  INTEGER  nelec

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Indexvariable
  INTEGER  i

!!!$.....................................................................

!!!$     Konstantenvektor auf Null setzen
!!$  b_komp = 0.
  do i=1,sanz
     b(i,1) = CMPLX(0d0)
  end do

!!!$     Aufbau des Konstanten- bzw. Stromvektors mit Skalierung
!!!$     ( A * x + b = 0 )
  b(enr(nelec),1) = CMPLX(-1_prec)

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

  if (lsink) b(nsink,1) = CMPLX(1_prec)

  return
end subroutine kompb

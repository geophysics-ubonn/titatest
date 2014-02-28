subroutine scalab(a_scal,b_scal,fak_scal)

!!!$     Unterprogramm skaliert 'a' und 'b' und liefert die Skalierungs-
!!!$     faktoren im Vektor 'fak'.

!!!$     ( Vgl. Subroutine 'SCALBNDN' in Schwarz (1991) )

!!!$     Andreas Kemna                                            11-Oct-1993
!!!$     Letzte Aenderung   07-Mar-2003

!!!$.....................................................................

  USE alloci
  USE femmod
  USE elemmod
  USE errmod

  IMPLICIT none


!!!$.....................................................................
!!$Gesamtsteifigkeitsmatrix
  COMPLEX (prec),DIMENSION(*):: a_scal
!!$ Berechnete Potentialwerte (bzw. Loesungsverktor)
  COMPLEX (prec),DIMENSION(*):: b_scal
!!$ Skalirerungsfaktor
  REAL (prec),DIMENSION(*)  :: fak_scal

!!!$     Hilfsvariablen
  INTEGER (KIND=4) ::    idi,i0
  INTEGER (KIND=4) ::     ja
  REAL(prec)  ::     dum

!!!$     Indexvariablen
  INTEGER (KIND=4) ::     i,j

!!!$.....................................................................

  do  i=1,sanz

     idi = i*(mb+1)
     dum = ABS(a_scal(idi))

     if (dum.le.0d0) then
        WRITE (fetxt,*)'scalab idi',idi,'i',i
        errnr = 27
        goto 1000
     end if

     a_scal(idi) = a_scal(idi) / CMPLX(dum)
     fak_scal(i) = 1d0 / SQRT(dum)
     b_scal(i)   = b_scal(i) * CMPLX(fak_scal(i))

     if (i.eq.1) CYCLE

     i0 = i*mb
     ja = max0(1,i-mb)

     do j=ja,i-1
        a_scal(i0+j) = a_scal(i0+j) * CMPLX(fak_scal(i)*fak_scal(j))
     END do

  END do

  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

1000 return

end subroutine scalab

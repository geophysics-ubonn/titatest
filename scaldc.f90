subroutine scaldc(a_scal,b_scal,fak_scal)

!!!$     Unterprogramm skaliert 'adc' und 'bdc' und liefert die Skalierungs-
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
  REAL(KIND(0D0)),DIMENSION(*)  ::  a_scal
  REAL(KIND(0D0)),DIMENSION(*)  ::  b_scal
  REAL(KIND(0D0)),DIMENSION(*)  ::  fak_scal

!!!$     Hilfsvariablen
  INTEGER (KIND=4) ::     idi,i0
  INTEGER (KIND=4) ::     ja
  REAL(KIND(0D0))  ::     dum

!!!$     Indexvariablen
  INTEGER (KIND=4) ::     i,j

!!!$.....................................................................

  do i=1,sanz

     idi = i*(mb+1)
     dum = a_scal(idi)

     if (dum.le.0d0) then
        WRITE (fetxt,'(a,I6,A,I6)')'scaldc',i,'idi',idi
        errnr = 27
        RETURN
     end if

     a_scal(idi) = 1d0
     fak_scal(i)   = 1d0 / dsqrt(dum)
     b_scal(i)   = b_scal(i) * fak_scal(i)

     if (i.eq.1) CYCLE

     i0 = i*mb
     ja = max0(1,i-mb)

     do j=ja,i-1
        a_scal(i0+j) = a_scal(i0+j)*fak_scal(i)*fak_scal(j)
     end do

  END DO

  errnr = 0

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

end subroutine scaldc

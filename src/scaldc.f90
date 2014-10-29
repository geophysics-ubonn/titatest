SUBROUTINE scaldc(a_scal,b_scal,fak_scal)

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

  IMPLICIT NONE


!!!$.....................................................................
  REAL(KIND(0D0)),DIMENSION((mb+1)*sanz)  ::  a_scal
  REAL(KIND(0D0)),DIMENSION(sanz)  ::  b_scal
  REAL(KIND(0D0)),DIMENSION(sanz)  ::  fak_scal

!!!$     Hilfsvariablen
  INTEGER (KIND=4) ::     idi,i0
  INTEGER (KIND=4) ::     ja
  REAL(KIND(0D0))  ::     dum

!!!$     Indexvariablen
  INTEGER (KIND=4) ::     i,j

!!!$.....................................................................

  DO i=1,sanz

     idi = i*(mb+1)
     dum = a_scal(idi)

     IF (dum.LE.0d0) THEN
        WRITE (fetxt,'(a,I6,A,I6)')'scaldc',i,'idi',idi
        errnr = 27
        RETURN
     END IF

     a_scal(idi) = 1d0
     fak_scal(i)   = 1d0 / dsqrt(dum)
     b_scal(i)   = b_scal(i) * fak_scal(i)

     IF (i == 1) CYCLE

     i0 = i*mb
     ja = max0(1,i-mb)

     DO j=ja,i-1
        a_scal(i0+j) = a_scal(i0+j)*fak_scal(i)*fak_scal(j)
     END DO

  END DO

  errnr = 0

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen
END SUBROUTINE scaldc

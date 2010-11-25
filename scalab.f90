subroutine scalab()

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

!!!$     Hilfsvariablen
  INTEGER (KIND=4) ::    idi,i0
  INTEGER (KIND=4) ::     ja
  REAL(KIND(0D0))  ::     dum

!!!$     Indexvariablen
  INTEGER (KIND=4) ::     i,j

!!!$.....................................................................

  do  i=1,sanz

     idi = i*(mb+1)
     dum = cdabs(a(idi))

     if (dum.le.0d0) then
        WRITE (fetxt,*)'scalab idi',idi,'i',i
        errnr = 27
        goto 1000
     end if

     a(idi) = a(idi) / dcmplx(dum)
     fak(i) = 1d0 / dsqrt(dum)
     b(i)   = b(i) * dcmplx(fak(i))

     if (i.eq.1) CYCLE

     i0 = i*mb
     ja = max0(1,i-mb)

     do j=ja,i-1
        a(i0+j) = a(i0+j) * dcmplx(fak(i)*fak(j))
     END do

  END do

  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

1000 return

end subroutine scalab

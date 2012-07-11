subroutine choldc(a_chol)

!!!$  Cholesky-Zerlegung der positiv definiten Matrix 'adc'; erfolgt auf dem
!!!$  Platz von 'adc', d.h. bei Auftreten eines Fehlers ist gegebene Matrix
!!!$    'adc' zerstoert.

!!!$   ( Vgl. Subroutine 'CHOBNDN' in Schwarz (1991) )

!!!$   Andreas Kemna                                            11-Oct-1993
!!!$   Letzte Aenderung   07-Mar-2003

!!!$.....................................................................

  USE alloci
  USE elemmod
  USE errmod

  IMPLICIT none


!!!$.....................................................................
  REAL (KIND(0D0)),DIMENSION(*)    ::  a_chol

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsvariablen
  INTEGER (KIND = 4)  ::     idi,i0,ij,j0
  INTEGER (KIND = 4)  ::     m1,fi
  REAL (KIND(0D0))    ::     s

!!!$     Indexvariablen
  INTEGER (KIND = 4)  ::     i,j,k

!!!$.....................................................................

  m1 = mb+1

  do i=1,sanz

     idi = i*m1
     fi  = max0(1,i-mb)
     i0  = idi-i

     do j=fi,i

        ij = i0+j
        j0 = j*mb
        s  = a_chol(ij)

        do k=fi,j-1
           s = s - a_chol(i0+k)*a_chol(j0+k)
        END do

        if (j.lt.i) then

           a_chol(ij) = s / a_chol(j*m1)

        else

           if (s.le.0d0) then
              fetxt = ' '
              errnr = 28
              goto 1000
           end if

           a_chol(idi) = dsqrt(s)

        end if

     END DO

  END do

  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

1000 return

end subroutine choldc
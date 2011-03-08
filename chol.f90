subroutine chol()

!!!$     Cholesky-Zerlegung der positiv definiten Matrix 'a'; erfolgt auf dem
!!!$     Platz von 'a', d.h. bei Auftreten eines Fehlers ist gegebene Matrix
!!!$     'a' zerstoert.

!!!$     ( Vgl. Subroutine 'CHOBNDN' in Schwarz (1991) )

!!!$     Andreas Kemna                                            11-Oct-1993
!!!$     Letzte Aenderung   07-Mar-2003

!!!$.....................................................................

  USE alloci
  USE elemmod
  USE errmod
  USE ompmod

  IMPLICIT none


!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsvariablen
  INTEGER (KIND = 4)  ::     idi,i0,ij,j0
  INTEGER (KIND = 4)  ::     m1,fi
  COMPLEX (KIND(0D0)) ::  s

!!!$     Indexvariablen
  INTEGER (KIND = 4)  ::     i,j,k

!!!$.....................................................................

  m1 = mb+1

  !$OMP PARALLEL PRIVATE(idi,i0,ij,j0,s,k,i,j,fi) 
  !$OMP DO ORDERED
  do i=1,sanz
     !$OMP ORDERED
     idi = i*m1
     fi  = max0(1,i-mb)
     i0  = idi-i
     
     do j=fi,i
        
        ij = i0+j
        j0 = j*mb
        s  = a(ij)

        do k=fi,j-1
           s = s - a(i0+k)*a(j0+k)
        END do

        if (j.lt.i) then

           a(ij) = s / a(j*m1)

        else

           if (cdabs(s).le.0d0) then
              fetxt = ' '
              errnr = 28
              STOP
           end if

           a(idi) = cdsqrt(s)

        end if

     END do
     !$OMP END ORDERED
  END do
  !$OMP END DO
  !$OMP END PARALLEL
  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

1000 return

end subroutine chol

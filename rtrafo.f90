subroutine rtrafo()

!!!$     Unterprogramm zur Ruecktransformation.

!!!$     Andreas Kemna                                            20-Dec-1993
!!!$     Letzte Aenderung   13-Nov-1997

!!!$.....................................................................

  USE alloci
  USE femmod
  USE electrmod
  USE elemmod
  USE wavenmod

  IMPLICIT none

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Pi
  REAL (KIND(0D0))  ::   pi

!!!$     Hilfsvariablen
  COMPLEX (KIND(0D0))  ::    summe
  REAL (KIND(0D0))    ::     summdc

!!!$     Indexvariablen
  INTEGER (KIND=4)    ::j,k,l

!!!$.....................................................................

  pi = dacos(-1d0)

  if (ldc) then

     do l=1,eanz
        do j=1,sanz
           summdc = 0d0

           do k=1,kwnanz
              summdc = summdc + kpotdc(j,l,k)*kwnwi(k)
           end do

           hpotdc(j,l) = summdc / pi
        end do
     end do

  else

     do l=1,eanz
        do j=1,sanz
           summe = dcmplx(0d0)

           do k=1,kwnanz
              summe = summe + kpot(j,l,k)*dcmplx(kwnwi(k))
           end do

           hpot(j,l) = summe / dcmplx(pi)
        end do
     end do

  end if

  return
end subroutine rtrafo

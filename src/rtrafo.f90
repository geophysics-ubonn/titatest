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
  REAL (prec)  ::   pi

!!!$     Hilfsvariablen
  COMPLEX (prec)  ::    summe
  REAL (prec)    ::     summdc

!!!$     Indexvariablen
  INTEGER (KIND=4)    ::j,k,l

!!!$.....................................................................
  print*,'inverse Fourier transform'
  pi = dacos(-1d0)

  if (ldc) then

     !$OMP PARALLEL DEFAULT (none) &
     !$OMP PRIVATE (summdc,l,j,k) &
     !$OMP SHARED (eanz,sanz,kwnanz,kwnwi,kpotdc,hpotdc,pi)
     !$OMP DO COLLAPSE (2) PRIVATE (l,j,k,summdc)
     do l=1,eanz
        do j=1,sanz
           summdc = 0d0

           do k=1,kwnanz
              summdc = summdc + kpotdc(j,l,k)*kwnwi(k)
           end do

           hpotdc(j,l) = summdc / pi
!!$           IF (summdc < EPSILON(0D0)) THEN
!!$              print*,'rtrafo::',l,j
!!$              STOP
!!$           END IF
        end do
     end do
     !$OMP END PARALLEL

  else

     !$OMP PARALLEL DEFAULT (none) &
     !$OMP PRIVATE(summe,l,j,k) &
     !$OMP SHARED (eanz,sanz,kwnanz,kwnwi,kpot,hpot,pi)
     !$OMP DO COLLAPSE (2) PRIVATE (l,j,k,summe)
     do l=1,eanz
        do j=1,sanz
           summe = CMPLX(0d0)

           do k=1,kwnanz
              summe = summe + kpot(j,l,k)*CMPLX(kwnwi(k))
           end do

           hpot(j,l) = summe / CMPLX(pi)
        end do
     end do
     !$OMP END PARALLEL

  end if

end subroutine rtrafo

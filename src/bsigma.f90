SUBROUTINE bsigma
!!!$     
!!!$     Unterprogramm zur Belegung des Modell vektors (sigma)
!!!$     mit dem verbesserten Modell aus par
!!!$     im Moment noch unspektakulaer, da elanz = manz 
!!!$     
!!!$     Copyright by Andreas Kemna
!!!$     
!!!$     Erste Version von Roland Martin                          19-Mar-2010
!!!$     
!!!$     Letzte Aenderung                                         23-Mar-2010
!!!$     
!!!$.....................................................................

  USE invmod   ! fuer par
  USE sigmamod ! fuer sigma
  USE modelmod ! fuer manz und mnr
  USE elemmod  ! fuer elanz
  USE errmod ! fuer fetxt und errnr

  IMPLICIT none

!!!$.....................................................................

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsfeld
!!!$      LOGICAL,DIMENSION(:),ALLOCATABLE :: lfeld
!!!$     Indexvariablen
  INTEGER (KIND = 4) :: i
!!!$.....................................................................


  errnr = 4

  do i=1,elanz
     sigma(i) = CDEXP(par(mnr(i)))
  end do

  errnr = 0

end SUBROUTINE bsigma

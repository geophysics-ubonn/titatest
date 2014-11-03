SUBROUTINE bpar
!!!$     
!!!$     Unterprogramm zur Belegung des Parameter Vektors
!!!$     im Moment noch unspektakulaer, da elanz = manz 
!!!$     
!!!$     Copyright by Andreas Kemna
!!!$     
!!!$     Erste Version von Roland Martin                          19-Mar-2010
!!!$     
!!!$     Letzte Aenderung                                         23-Mar-2010
!!!$     
!!!$.....................................................................

  USE invmod, ONLY: par   ! fuer par
  USE sigmamod,ONLY: sigma ! fuer sigma
  USE modelmod,ONLY:mnr ! fuer manz und mnr
  USE elemmod,ONLY:elanz  ! fuer elanz
  USE errmod,ONLY:errnr ! fuer fetxt und errnr

  IMPLICIT none

!!!$.....................................................................

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsfeld
!!!$      LOGICAL,DIMENSION(:),ALLOCATABLE :: lfeld
!!!$     Indexvariablen
  INTEGER (KIND = 4) ::     i
!!!$.....................................................................


  errnr = 4

  do i=1,elanz
     par(mnr(i))   = CDLOG(sigma(i))
  end do

  errnr = 0

end SUBROUTINE bpar

      SUBROUTINE bresdc
c     
c     Unterprogramm berechnet Aufloesungsmatrix
c     (A^TC_d^-1A + C_m^-1)^-1 RES = A^TC_d^-1A
c     Fuer beliebige Triangulierung
c     
c     Andreas Kemna                                            02-Nov-2009
c     
c     Letzte Aenderung    RM                                   06-Nov-2009
c     
c.........................................................................
      USE alloci
      
      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'konv.fin'
      INCLUDE 'model.fin'
      INCLUDE 'err.fin'
!.....................................................................
!     PROGRAMMINTERNE PARAMETER:
!     Hilfsvariablen 
      INTEGER                                      :: i
      REAL(KIND(0D0)),DIMENSION(:,:),ALLOCATABLE   :: work
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE     :: ipiv
!.....................................................................

c$$$  solve (A^TC_d^-1A + C_m^-1) x = A^TC_d^-1A

      ALLOCATE (work(manz,manz))
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem work in bresdc ',
     1        manz**2.*8/(1024**3.),' GB'
         errnr = 97
         RETURN
      END IF
      ALLOCATE (ipiv(manz))
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem IPIV in bresdc'
         errnr = 97
         RETURN
      END IF

      work = atadcreg ! work is replaced by PLU decomposition

c$$$  setting up RHS, overwriting atadcreg
      atadcreg = atadc

c$$$  Solving Linear System Ax=B -> x=A^-1B

      WRITE (*,'(a)')ACHAR(9)//'Solving Ax=B'
      CALL DGESV(manz,manz,work,manz,ipiv,atadcreg,manz,errnr)
      
      IF (errnr /= 0) THEN
         PRINT*,'Zeile::',atadcreg(abs(errnr),:)
         PRINT*,'Spalte::',atadcreg(:,abs(errnr))
         errnr = 108
         RETURN
      END IF
      
      DEALLOCATE (work,ipiv)

      END

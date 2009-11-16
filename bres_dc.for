      SUBROUTINE bres_dc(kanal,ols)
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
      INTEGER                                      :: i,kanal
      REAL(KIND(0D0)),DIMENSION(:,:),ALLOCATABLE   :: work
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE     :: ipiv
      LOGICAL,OPTIONAL                             :: ols 
!     switch solv ordinary linear system or not^^
!.....................................................................

c$$$  solve (A^TC_d^-1A + C_m^-1) x = A^TC_d^-1A

      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4

      IF (.NOT.PRESENT(ols)) THEN
         ata_reg_dc = MATMUL(cov_m_dc,ata_dc)
      ELSE
         ALLOCATE (work(manz,manz),STAT=errnr)
         IF (errnr/=0) THEN
            WRITE (*,'(/a,G10.3,a/)')'Allocation problem work'//
     1           ' in bres_dc',REAL(manz)**2.*8./(1024.**3.),' GB'
            errnr = 97
            RETURN
         END IF
         ALLOCATE (ipiv(manz),STAT=errnr)
         IF (errnr/=0) THEN
            WRITE (*,'(/a/)')'Allocation problem IPIV in bresdc'
            errnr = 97
            RETURN
         END IF
         work = ata_reg_dc      ! work is replaced by LU decomposition
c$$$  setting up RHS, overwriting atadcreg
         ata_reg_dc = ata_dc
c$$$  Solving Linear System Ax=B -> x=A^-1B
         WRITE (*,'(a)')ACHAR(9)//'Solving Ax=B'
         CALL DGESV(manz,manz,work,manz,ipiv,ata_reg_dc,manz,errnr)
         IF (errnr /= 0) THEN
            PRINT*,'Zeile::',ata_reg_dc(abs(errnr),:)
            PRINT*,'Spalte::',ata_reg_dc(:,abs(errnr))
            errnr = 108
            RETURN
         END IF
         DEALLOCATE (work,ipiv)
      END IF

      WRITE (kanal,*)manz
      DO i=1,manz
         WRITE (kanal,*)log10(abs(ata_reg_dc(i,i))),ata_reg_dc(i,i)
      END DO

      CLOSE (kanal)
      errnr = 0
 999  RETURN

      END

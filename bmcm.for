      SUBROUTINE bmcm(kanal)
c     
c     Unterprogramm berechnet (einfache) Modell Kovarianz Matrix
c     A^TC_d^-1A + C_m^-1
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
      INCLUDE 'elem.fin'
      INCLUDE 'err.fin'
!.....................................................................
!     PROGRAMMINTERNE PARAMETER:
!     Hilfsvariablen 
      INTEGER                                       :: i,kanal
      COMPLEX(KIND(0D0)),DIMENSION(:,:),ALLOCATABLE :: work
      COMPLEX(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: ipiv
!.....................................................................

c$$$  solve A^TC_d^-1A + C_m^-1
      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4

      ALLOCATE (work(manz,manz),STAT=errnr)
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem WORK in bmcm'
         errnr = 97
         RETURN
      END IF
      ALLOCATE (ipiv(manz),STAT=errnr)
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem IPIV in bmcmdc'
         errnr = 97
         RETURN
      END IF
      
      work = ata_reg ! work is replaced by PLU decomposition
      
c$$$  building Right Hand Side (unit matrix)
      DO i=1,manz
         cov_m(i,i) = CMPLX(1.d0,1.d0)
      END DO

c$$$  Solving Linear System Ax=B -> B=A^-1
      WRITE (*,'(a)')'Solving Ax=B'
      CALL ZGESV(manz,manz,work,manz,ipiv,cov_m,manz,errnr)

      IF (errnr /= 0) THEN
         PRINT*,'Zeile::',cov_m(abs(errnr),:)
         PRINT*,'Spalte::',cov_m(:,abs(errnr))
         errnr = 108
         RETURN
      END IF

      DEALLOCATE (work,ipiv)
      WRITE (kanal,*)manz
      DO i=1,manz
         WRITE (kanal,*)log10(sqrt(abs(REAL(cov_m(i,i))))),
     1        log10(sqrt(abs(dimag(cov_m(i,i)))))
      END DO

 999  RETURN
      END

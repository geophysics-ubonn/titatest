      SUBROUTINE bmcm
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
      INTEGER                                       :: i,j,k
      COMPLEX(KIND(0D0)),DIMENSION(:,:),ALLOCATABLE :: work
      COMPLEX(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: ipiv
!.....................................................................

c$$$  solve A^TC_d^-1A + C_m^-1

      ALLOCATE (ata(manz,manz),STAT=errnr)
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem ATA in bmcm'
         errnr = 97
         RETURN
      END IF
      ALLOCATE (work(manz,manz))
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem WORK in bmcm'
         errnr = 97
         RETURN
      END IF
      ALLOCATE (cov_m(manz,manz))
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem MCM_1 in bmcm'
         errnr = 97
         RETURN
      END IF
      ALLOCATE (ipiv(manz))
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem IPIV in bmcmdc'
         errnr = 97
         RETURN
      END IF

      ata = MATMUL(TRANSPOSE(sens),sens) ! is stored for later
      
      work = ata              ! preselect A^TC_d^-1A 
      IF (ltri == 0) THEN
         DO j=1,manz
            DO i=1,manz
               IF (i==j) THEN   ! sparse C_m
                  work(i,j) = work(i,j) + lam *
     1                 CMPLX(smatm(i,1),smatm(i,1))
                  IF (i+1 < manz) work(i+1,j) = work(i+1,j) + lam * 
     1                 CMPLX(smatm(i+1,2),smatm(i+1,2))
                  IF (i+nx < manz) work(i+1,j) = work(i+1,j) + lam *
     1                 CMPLX(smatm(i+nx,3),smatm(i+nx,3))
                  IF (i-1 > 1) work(i-1,j) = work(i-1,j) + lam *
     1                 CMPLX(smatm(i-1,2),smatm(i-1,2))
                  IF (i-nx > 1) work(i-nx,j) = work(i-nx,j) + lam *
     1                 CMPLX(smatm(i-nx,3),smatm(i-nx,3))
               END IF
            END DO
         END DO
      ELSE IF (ltri < 10 ) THEN
         DO j=1,manz
            DO i=1,manz
               IF (i==j) THEN   ! sparse C_m
                  work(i,j) = work(i,j) + lam *
     1                 CMPLX(smatm(i,0),smatm(i,0))
                  DO k=1,nachbar(i,0)
                     IF (nachbar(i,k) /= 0) work(nachbar(i,k),j) = 
     1                    work(nachbar(i,k),j) + lam *
     1                    CMPLX(smatm(i,k),smatm(i,k))
                  END DO
               END IF
            END DO
         END DO
      ELSE
         work = work + CMPLX(smatm,smatm) ! for full C_m..
      END IF

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

      END

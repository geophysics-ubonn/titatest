      subroutine bsmatmtv      !betatv
c     
c     Unterprogramm belegt die Rauhigkeitsmatrix mit total variance
c     fuer beliebige Triangulierung 
c
c     Copyright by Andreas Kemna 2009
c     
c     Erste Version von Roland Martin                          23-Nov-2009
c     
c     Letzte Aenderung   RM                                    23-Nov-2009
c     
c.........................................................................
      USE alloci
      USE femmod
      USE datmod
      
      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'model.fin'       ! mit nachbar und ldir
      INCLUDE 'konv.fin'
      INCLUDE 'inv.fin'
!.....................................................................
!     PROGRAMMINTERNE PARAMETER:
!     Hilfsvariablen 
      REAL(KIND(0D0)) :: dum
      INTEGER         :: i,j,l,k,smaxs,ik,anz
      REAL(KIND(0D0)) :: edglen ! Kantenlaenge
      REAL(KIND(0D0)) :: dist   ! Abstand der Schwerpunkte
      REAL(KIND(0D0)) :: sp1(2),sp2(2) ! Schwerpunktkoordinaten
      REAL(KIND(0D0)) :: ang    !Winkel fuer anisotrope Glaettung
      REAL(KIND(0D0)) :: snsmn  !Mittlere Sensitivität
      REAL(KIND(0D0)) :: alfgeo !Anisotrope (geometrische) Glaettung
      REAL(KIND(0D0)) :: alftv  !TV Glaettung
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: csens 
!.....................................................................
      
      IF (.NOT.ALLOCATED(csens)) ALLOCATE (csens(manz))
      csens=0.

      IF (lip) THEN
         DO j=1,manz
            DO i=1,nanz
               csens(j) = csens(j) + 
     1              DBLE(sens(i,j))*DBLE(sens(i,j))* 
     1              wmatd(i)*DBLE(wdfak(i))
            END DO
         END DO
      ELSE IF (ldc) THEN
         DO j=1,manz
            DO i=1,nanz
               csens(j) = csens(j) + 
     1              sensdc(i,j)*sensdc(i,j)* 
     1              wmatd(i)*DBLE(wdfak(i))
            END DO
         END DO
      ELSE
         DO j=1,manz
            DO i=1,nanz
               csens(j) = csens(j) + 
     1              DCONJG(sens(i,j))*sens(i,j)* 
     1              wmatd(i)*dble(wdfak(i)) ! wechselt automatisch zu
!     wmatdp bei lip
            END DO
         END DO
      ENDIF
      dum = MAXVAL(csens)
      IF (dum.gt.1d-12) THEN
         csens = csens/dum
      END IF
      snsmn = 0.
      DO i=1,manz
         snsmn = snsmn + csens(i)
      END DO
      snsmn = snsmn / DBLE(manz)
!     Summe der Sensitivitaeten normieren
      WRITE(*,*) 'dum snsmn',dum,snsmn

      smaxs=MAXVAL(selanz)
      IF (.NOT.ALLOCATED(smatm)) ALLOCATE (smatm(manz,smaxs+1))
      smatm = 0d0               ! initialize smatm

      DO i=1,elanz

         sp1 = 0.
         
         DO k=1,smaxs           ! Schwerpunkt berechnen 
            sp1(1) = sp1(1) + sx(snr(nrel(i,k)))
            sp1(2) = sp1(2) + sy(snr(nrel(i,k)))
         END DO
         
         sp1(1) = sp1(1)/smaxs
         sp1(2) = sp1(2)/smaxs ! Mittelpunkt des aktuellen Elements
         
         DO k=1,smaxs           ! jedes flaechenele hat mind einen nachbarn

            ik = MOD(k,smaxs) + 1

            edglen = SQRT((sx(snr(nrel(i,k))) - 
     1           sx(snr(nrel(i,ik))))**2 +
     1           (sy(snr(nrel(i,k))) -
     1           sy(snr(nrel(i,ik))))**2) ! edge
            

            IF (nachbar(i,k)>0) THEN !nachbar existiert 
               
               sp2 = 0.

               DO l=1,smaxs
                  sp2(1) = sp2(1) + sx(snr(nrel(nachbar(i,k),l)))
                  sp2(2) = sp2(2) + sy(snr(nrel(nachbar(i,k),l)))
               END DO
               
               sp2(1) = sp2(1)/smaxs ! schwerpunkt des nachbar elements
               sp2(2) = sp2(2)/smaxs
               
!     Geometrischer Teil...
               dist = SQRT((sp1(1) - sp2(1))**2. + 
     1              (sp1(2) - sp2(2))**2.)

               ang = DATAN2((sp1(2) - sp2(2)),(sp1(1) - sp2(1))) !neu

               alfgeo = DSQRT((alfx*DCOS(ang))**2. + 
     1              (alfz*DSIN(ang))**2.)

               alftv = edglen / dist * alfgeo
!     Total variance
               dum = SQRT(alftv**2. + betamgs**2.)
!     nun glaettung belegen
               smatm(i,k) = -dum ! neben Diagonale
               smatm(i,smaxs+1) = smatm(i,smaxs+1) + dum !Hauptdiagonale
            END IF

         END DO
      END DO

      END


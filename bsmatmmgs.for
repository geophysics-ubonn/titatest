      subroutine bsmatmmgs      !MGS
c     
c     Unterprogramm belegt die Rauhigkeitsmatrix ala Portniaguine und Zhdanov [1999]
c     Fuer beliebige Triangulierung mit Sensitivitäten gewichtet [Blaschek 2008]
c
c     Copyright by Andreas Kemna                               2009
c     
c     Created by Roland Blaschek/Roland Martin                 03-Nov-2009
c     
c     Last edited  RM                                          18-Dec-2009
c     
c.........................................................................
      USE alloci
      
      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'dat.fin'
      INCLUDE 'model.fin'       ! mit nachbar und ldir
      INCLUDE 'konv.fin'
      INCLUDE 'inv.fin'
      INCLUDE 'fem.fin'
!.....................................................................
!     PROGRAMMINTERNE PARAMETER:
!     Hilfsvariablen 
      REAL(KIND(0D0)) :: dum
      INTEGER         :: i,j,l,k,smaxs,ik,anz
      REAL(KIND(0D0)) :: edglen(selmax) ! Kantenlaenge
      REAL(KIND(0D0)) :: dist(selmax) ! Abstand der Schwerpunkte
      REAL(KIND(0D0)) :: sp(0:selmax,2) ! Schwerpunktkoordinaten
      REAL(KIND(0D0)) :: ang    !Winkel fuer anisotrope Glaettung
      REAL(KIND(0D0)) :: snsmn  !Mittlere Sensitivität
      REAL(KIND(0D0)) :: alfgeo !Anisotrope Glaettung
      REAL(KIND(0D0)) :: alfmgs !MGS Glaettung
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: csens 
!.....................................................................
      
      smaxs=MAXVAL(selanz)

      IF (.NOT.ALLOCATED(smatm)) 
     1     ALLOCATE (smatm(manz,smaxs+1),STAT=errnr)

      IF (.NOT.ALLOCATED(csens)) 
     1     ALLOCATE (csens(manz),STAT=errnr)

      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem bsmatmmgs'
         
         errnr = 97
         RETURN
      END IF

      csens=0.
      smatm = 0d0               ! initializing

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
 
      DO i=1,elanz

         sp(0:smaxs,:) = 0.
         
         DO k=1,smaxs           ! Schwerpunkt berechnen 
            sp(0,1) = sp(0,1) + sx(snr(nrel(i,k)))
            sp(0,2) = sp(0,2) + sy(snr(nrel(i,k)))
         END DO
         
         sp(0,1) = sp(0,1)/smaxs
         sp(0,2) = sp(0,2)/smaxs ! Mittelpunkt des aktuellen Elements
         
         DO k=1,smaxs           ! jedes flaechenele hat mind einen nachbarn

            ik = MOD(k,smaxs) + 1

            edglen(k) = SQRT((sx(snr(nrel(i,k))) - 
     1           sx(snr(nrel(i,ik))))**2 +
     1           (sy(snr(nrel(i,k))) -
     1           sy(snr(nrel(i,ik))))**2) ! edge
            

            IF (nachbar(i,k)>0) THEN !nachbar existiert 
               
               DO l=1,smaxs
                  sp(k,1) = sp(k,1) + sx(snr(nrel(nachbar(i,k),l)))
                  sp(k,2) = sp(k,2) + sy(snr(nrel(nachbar(i,k),l)))
               END DO
               
               sp(k,1) = sp(k,1)/smaxs ! schwerpunkt des nachbar elements
               sp(k,2) = sp(k,2)/smaxs
               
!     Geometrischer Teil...
               dist(k) = SQRT((sp(0,1) - sp(k,1))**2 +
     1              (sp(0,2) - sp(k,2))**2)
               ang = DATAN2((sp(0,2) - sp(k,2)),(sp(0,1) - sp(k,1))) !neu
               alfgeo = DSQRT((alfx*DCOS(ang))**2. + 
     1              (alfz*DSIN(ang))**2.)
!     MGS Teil
               dum = CDABS(par(i)-par(nachbar(i,k)))
               IF (ltri == 2) THEN ! reines MGS
                  dum = dum 
               ELSE IF (ltri == 3) THEN ! sensitivitaetswichtung 1
                  dum = dum * (1d0 + 0.2d0 * (DABS( DLOG10(csens(i)) + 
     1                 DLOG10(csens(nachbar(i,k))) ) ))
               ELSE IF (ltri == 4) THEN ! sensitivitaetswichtung und mittelwert
                  dum = dum * (1d0 + 0.2d0 * (DABS( DLOG10(csens(i)) + 
     1                 DLOG10(csens(nachbar(i,k))) ) / snsmn ))
               END IF
               alfmgs = 1d0 - dum**2. / (dum**2. + betamgs**2.)
!     gesamt eintrag
               dum = edglen(k) / dist(k) * alfgeo * alfmgs
!     nun glaettung belegen
               smatm(i,k) = -dum ! neben Diagonale
               smatm(i,smaxs+1) = smatm(i,smaxs+1) + dum !Hauptdiagonale
            END IF

         END DO
      END DO
      errnr = 0
      END


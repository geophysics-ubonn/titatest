      subroutine bsmatmmgs      !MGS
c     
c     Unterprogramm belegt die Rauhigkeitsmatrix ala Portniaguine und Zhdanov [1999]
c     Fuer beliebige Triangulierung mit Sensitivitäten gewichtet [Blaschek 2008]
c
c     Copyright by Andreas Kemna 2009
c     
c     Erste Version von Roland Martin                          03-Nov-2009
c     
c     Last edited  RM                                          18-Dec-2009
c     
c.........................................................................

      USE alloci
      USE femmod
      USE datmod
      USE invmod
      USE sigmamod       ! for sigma
      
      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'model.fin'       ! mit nachbar und ldir
      INCLUDE 'konv.fin'
      INCLUDE 'err.fin'

!.....................................................................
!     PROGRAMMINTERNE PARAMETER:
!     Hilfsvariablen 
      REAL(KIND(0D0)) :: dum,dum2
      INTEGER         :: i,j,l,k,smaxs,ik,anz
      REAL(KIND(0D0)) :: edglen ! Kantenlaenge
      REAL(KIND(0D0)) :: dist ! Abstand der Schwerpunkte
      REAL(KIND(0D0)) :: sp1(2),sp2(2) ! Schwerpunktkoordinaten
      REAL(KIND(0D0)) :: ang    !Winkel fuer anisotrope Glaettung
      REAL(KIND(0D0)) :: snsmn  !Mittlere Sensitivität
      REAL(KIND(0D0)) :: alfgeo !Anisotrope Glaettung
      REAL(KIND(0D0)) :: alfmgs !MGS Glaettung
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: csens 
!.....................................................................
      
      errnr = 4

      IF (.NOT.ALLOCATED(csens)) ALLOCATE (csens(manz),STAT=errnr)
      IF (errnr/=0) THEN
         fetxt = 'Allocation problem csens in bsmatmmgs'
         WRITE (*,'(/a/)')TRIM(fetxt)
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

! auf eins normieren...
      dum = MAXVAL(csens)
      IF (dum.gt.1d-12) csens = csens/dum

!     evtl auf Summer der Sensitivitaeten normieren
      snsmn = SUM (csens) / DBLE(manz)

c$$$      DO i=1,manz
c$$$         snsmn = snsmn + csens(i)
c$$$      END DO
c$$$      snsmn = snsmn / DBLE(manz)

      WRITE(*,*) 'dum snsmn',dum,snsmn

      smaxs=MAXVAL(selanz) ! triangles or rectangles
      IF (.NOT.ALLOCATED(smatm)) ALLOCATE (smatm(manz,smaxs+1),
     1     STAT=errnr)

      IF (errnr/=0) THEN
         fetxt = 'Allocation problem WORK in bmcm'
         WRITE (*,'(/a/)')TRIM(fetxt)
         errnr = 97
         RETURN
      END IF

      smatm = 0d0               ! initialize smatm

      DO i=1,elanz ! elanz = flaecheneles

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
     1           sy(snr(nrel(i,ik))))**2) ! edge of i,k and the next..
            

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
!     MGS Teil
               dum = CDABS(sigma(i) - sigma(nachbar(i,k))) / dist
!     Modell nur im zaehler nicht im Nenner -> fred fragen

!     \int \frac{(\nabla m_{ij})^2}{(\nabla m_{ij})^2+\beta^2}\;dA
!     -> (m_i-m_{i+1})^2 \frac{\Delta z_i}{\Delta x_i}
!     wobei (m_i-m_{i+1})^2 rausgezogen wird und spaeter 
!     als Matrix Vektor Produkt berechnet wird
!     (m_i-m_{i+1})^2 \frac{\Delta z_i}{\Delta x_i} 
!     -> smatm(i) = \frac{\Delta z_i}{\Delta x_i} * geometrischem Teil 
!     von anisotroper Wichtung
               IF (ltri == 5) THEN ! reines MGS

                  dum = dum**2. + betamgs**2.
                  dum = alfgeo * edglen / dist / dum

               ELSE IF (ltri == 6) THEN ! sensitivitaetswichtung 1 von RM
!     f(i,k) = 1 + g(i) + g(k)
                  dum2 = 1d0 + DABS(DLOG10(csens(i))) + 
     1                 DABS(DLOG10(csens(nachbar(i,k))))
!     dum2 = f(i,k)^2
                  dum2 = dum2**2.
!     dum = grad(m)^2 + (\beta/f(i,k)^2)^2
                  dum = dum**2. + (betamgs / dum2)**2.
!     dum = \alpha_{xz} * \Delta z / \Delta x / f(i,k)^2 / 
!     grad(m)^2 + (\beta/f(i,k)^2)^2
                  dum = alfgeo * edglen / dist /dum2 / dum

               ELSE IF (ltri == 7) THEN ! sensitivitaetswichtung 1 von RM
                  
!     f(i,k) = 1 + (g(i) + g(k))/mean(g)
                  dum2 = 1d0 + DABS((DLOG10(csens(i))) + 
     1                 DABS(DLOG10(csens(nachbar(i,k))))) / snsmn
!     dum2 = f(i,k)^2
                  dum2 = dum2**2.
!     dum = grad(m)^2 + (\beta/f(i,k)^2)^2
                  dum = dum**2. + (betamgs / dum2)**2.
!     dum = \alpha_{xz} * \Delta z / \Delta x / f(i,k)^2 / 
!     grad(m)^2 + (\beta/f(i,k)^2)^2
                  dum = alfgeo * edglen / dist /dum2 / dum

               ELSE IF (ltri == 8) THEN ! sensitivitaetswichtung von RB

!     der folgende code wurde mir so ueberliefert... 
!     keine ahnung was das genau macht
                  dum = dum * (1d0 + 0.2d0 * (DABS( DLOG10(csens(i)) + 
     1                 DLOG10(csens(nachbar(i,k))) ) ))
                  alfmgs = 1d0 - dum**2. / (dum**2. + betamgs**2.)
                  dum =  edglen * alfgeo * alfmgs
                  
               ELSE IF (ltri == 9) THEN
                  
                  dum = dum * (1d0 + 0.2d0 * (DABS( DLOG10(csens(i)) + 
     1                 DLOG10(csens(nachbar(i,k))) ) / snsmn ))
                  alfmgs = 1d0 - dum**2. / (dum**2. + betamgs**2.)
                  dum =  edglen * alfgeo * alfmgs

               END IF

!     nun glaettung belegen
               smatm(i,k) = -dum ! neben Diagonale
               smatm(i,smaxs+1) = smatm(i,smaxs+1) + dum !Hauptdiagonale

            END IF
         END DO
      END DO

      DEALLOCATE (csens)

      errnr = 0
 999  RETURN

      END


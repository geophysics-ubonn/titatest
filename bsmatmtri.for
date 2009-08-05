      subroutine bsmatmtri      !tri
c
c     Unterprogramm belegt die Rauhigkeitsmatrix....
c     Fuer beliebige Triangulierung
c
c     Andreas Kemna                                            29-Feb-1996
c
c     Letzte Aenderung                                         29-Jul-2009
c
c.........................................................................

      USE alloci
      
      INCLUDE 'parmax.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'dat.fin'
      INCLUDE 'model.fin'       ! mit nachbar und ldir
      INCLUDE 'konv.fin'
      INCLUDE 'inv.fin'
!.....................................................................

!     PROGRAMMINTERNE PARAMETER:

!     Hilfsvariablen 
      real            * 8     alfdis,dum
      integer         * 4     i,j,l,iflnr,ijdum,smaxs
      REAL            * 8     edglen(selmax)   ! Kantenlaenge
      REAL            * 8     dist(selmax)     ! Abstand der Schwerpunkte
      REAL            * 8     sp(0:selmax,2)   ! Schwerpunktkoordinaten
      REAL            * 8 ang   !Winkel fuer anisotrope Glaettung
      REAL            * 8 alfxz !Anisotrope Glaettung
      
!.....................................................................
      
      smaxs=MAXVAL(selanz)
      IF (.NOT.ALLOCATED(smatm)) ALLOCATE (smatm(manz,smaxs+1))
      smatm = 0d0               ! initialize smatm
      iflnr =0
      
      DO i=1,manz
         DO k=1,smaxs           ! Schwerpunkt des mittleren Elementes
            sp(0,1) = sp(0,1) + sx(snr(nrel(i,k)))
            sp(0,2) = sp(0,2) + sy(snr(nrel(i,k)))
         END DO
         
         sp(0,1) = sp(0,1)/smaxs
         sp(0,2) = sp(0,2)/smaxs ! Mittelpunkt des aktuellen Elements
         
         DO k=1,smaxs ! jedes flaechenele hat mind einen nachbarn

            ik = MOD(k,smaxs) + 1

            edglen(k) = SQRT((sx(snr(nrel(i,k))) - 
     1           sx(snr(nrel(i,ik))))**2 +
     1           (sy(snr(nrel(i,k))) -
     1           sy(snr(nrel(i,ik))))**2) ! edge
            
            sp(1:smaxs,1:2) = 0. ! mittlepunkte der Nachbarn initialisieren

            IF (nachbar(i,k)>0) THEN !nachbar existiert 
               
               DO l=1,smaxs
                  sp(k,1) = sp(k,1) + sx(snr(nrel(nachbar(i,k),l)))
                  sp(k,2) = sp(k,2) + sy(snr(nrel(nachbar(i,k),l)))
               END DO
               
               sp(k,1) = sp(k,1)/smaxs ! schwerpunkt des nachbar elements
               sp(k,2) = sp(k,2)/smaxs
               
               dist(k) = SQRT((sp(0,1) - sp(k,1))*(sp(0,1) - sp(k,1)) + 
     1              (sp(0,2) - sp(k,2))*(sp(0,2) - sp(k,2)))
               
               ang = DATAN2(sp(0,2) - sp(k,2),sp(0,1) - sp(k,1)) !neu
               
               alfxz = DSQRT((alfx*DCOS(ang))**2 + (alfz*DSIN(ang))**2)
               
               smatm(i,k) = -edglen(k)/dist(k)*alfxz
               
               smatm(i,smaxs+1) = smatm(i,smaxs+1) + 
     1              edglen(k)/dist(k)*alfxz !Hauptdiagonale
            END IF

         END DO
      END DO

c$$$      DO i=1,typanz
c$$$         
c$$$         IF (typ(i)>10) EXIT ! die Flächen kommen zuerst im gitter..
c$$$
c$$$         DO j=1,nelanz(i)
c$$$            iflnr = iflnr + 1 
c$$$            sp(0,1)=0. ; sp(0,2)=0.
c$$$            
c$$$            DO k=1,selanz(i)    ! Schwerpunkt des mittleren Elementes
c$$$               sp(0,1)=sp(0,1)+sx(snr(nrel(iflnr,k)))
c$$$               sp(0,2)=sp(0,2)+sy(snr(nrel(iflnr,k)))
c$$$            END DO
c$$$            
c$$$            sp(0,1)=sp(0,1)/selanz(i); sp(0,2)=sp(0,2)/selanz(i) ! Mittelpunkt des aktuellen Elements
c$$$            DO k=1,selanz(i)
c$$$
c$$$               edglen(k)=SQRT((sx(snr(nrel(iflnr,k)))-  ! gemeinsame Kante
c$$$     1              sx(snr(nrel(iflnr,MOD(k,selanz(i))+1))))**2+
c$$$     1              (sy(snr(nrel(iflnr,k)))- 
c$$$     1              sy(snr(nrel(iflnr,MOD(k,selanz(i))+1))))**2)
c$$$               
c$$$               sp(1:selmax,1:2)=0. ! mittlepunkte der Nachbarn initialisieren
c$$$
c$$$               IF (nachbar(iflnr,k)>0) THEN !nachbar existiert 
c$$$
c$$$                  ijdum=nachbar(nachbar(iflnr,k),0)
c$$$
c$$$                  DO l=1,ijdum  
c$$$
c$$$                     sp(k,1)=sp(k,1)+
c$$$     1                    sx(snr(nrel(nachbar(iflnr,k),l)))
c$$$                     sp(k,2)=sp(k,2)+
c$$$     1                    sy(snr(nrel(nachbar(iflnr,k),l)))
c$$$                  END DO
c$$$
c$$$                  sp(k,1)=sp(k,1)/ijdum ! schwerpunkt des nachbar elements
c$$$                  sp(k,2)=sp(k,2)/ijdum
c$$$
c$$$                  dist(k)=SQRT((sp(0,1)-sp(k,1))*(sp(0,1)-sp(k,1))+ 
c$$$     1                 (sp(0,2)-sp(k,2)) * (sp(0,2)-sp(k,2)))
c$$$
c$$$                  ang=DATAN2(sp(0,2)-sp(k,2),sp(0,1)-sp(k,1)) !neu
c$$$
c$$$                  alfxz=DSQRT((alfx*DCOS(ang))**2+
c$$$     1                 (alfz*DSIN(ang))**2)
c$$$
c$$$                  smatm(iflnr,k)=-edglen(k)/dist(k)*alfxz
c$$$                  
c$$$                  smatm(iflnr,selanz(1)+1)=smatm(iflnr,smaxs+1)+
c$$$     1                 edglen(k)/dist(k)*alfxz !0-Mitte
c$$$
c$$$               ENDIF
c$$$
c$$$            END DO 
c$$$         END DO 
c$$$      END DO
      
      END


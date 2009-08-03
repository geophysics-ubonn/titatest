      subroutine bsmatmtri      !tri
      
c     Unterprogramm belegt die Rauhigkeitsmatrix.
c
c     Andreas Kemna/Roland Blaschek                            29-Feb-1996
c     Letzte Aenderung                                         29-Jul-2009
c.........................................................................

      USE alloci
      
      INCLUDE 'parmax.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'dat.fin'
      INCLUDE 'model.fin'       ! mit nachbar und ldir
      INCLUDE 'konv.fin'
      INCLUDE 'inv.fin'
!     smatm(1:anztri,0:selmax) 
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
      DO i=1,typanz
         
         IF (typ(i)>10) EXIT ! die Flächen kommen zuerst im gitter..

         DO j=1,nelanz(i)
            iflnr = iflnr + 1 
            sp(0,1)=0. ; sp(0,2)=0.
            
            DO k=1,selanz(i)    ! Schwerpunkt des mittleren Elementes
               sp(0,1)=sp(0,1)+sx(snr(nrel(iflnr,k)))
               sp(0,2)=sp(0,2)+sy(snr(nrel(iflnr,k)))
            END DO
            
            sp(0,1)=sp(0,1)/selanz(i); sp(0,2)=sp(0,2)/selanz(i) ! Mittelpunkt des aktuellen Elements
            DO k=1,selanz(i)

               edglen(k)=SQRT((sx(snr(nrel(iflnr,k)))-  ! gemeinsame Kante
     1              sx(snr(nrel(iflnr,MOD(k,selanz(i))+1))))**2+
     1              (sy(snr(nrel(iflnr,k)))- 
     1              sy(snr(nrel(iflnr,MOD(k,selanz(i))+1))))**2)
               
               sp(1:selmax,1:2)=0. ! mittlepunkte der Nachbarn initialisieren

               IF (nachbar(iflnr,k)>0) THEN !nachbar existiert 

                  ijdum=nachbar(nachbar(iflnr,k),0)

                  DO l=1,ijdum  

                     sp(k,1)=sp(k,1)+
     1                    sx(snr(nrel(nachbar(iflnr,k),l)))
                     sp(k,2)=sp(k,2)+
     1                    sy(snr(nrel(nachbar(iflnr,k),l)))
                  END DO

                  sp(k,1)=sp(k,1)/ijdum ! schwerpunkt des nachbar elements
                  sp(k,2)=sp(k,2)/ijdum

                  dist(k)=SQRT((sp(0,1)-sp(k,1))*(sp(0,1)-sp(k,1))+ 
     1                 (sp(0,2)-sp(k,2)) * (sp(0,2)-sp(k,2)))

                  ang=DATAN2(sp(0,2)-sp(k,2),sp(0,1)-sp(k,1)) !neu

                  alfxz=DSQRT((alfx*DCOS(ang))**2+
     1                 (alfz*DSIN(ang))**2)

                  smatm(iflnr,k)=-edglen(k)/dist(k)*alfxz
                  
                  smatm(iflnr,selanz(1)+1)=smatm(iflnr,smaxs+1)+
     1                 edglen(k)/dist(k)*alfxz !0-Mitte

               ENDIF

            END DO 
         END DO 
      END DO
      
      END


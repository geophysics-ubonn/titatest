      subroutine bsmatm3()      !tri
      
! Unterprogramm belegt die Rauhigkeitsmatrix.

! Andreas Kemna                                            29-Feb-1996
!                                       Letzte Aenderung   23-Apr-1998
! Aenderung auf allgemeine Netze, Roland Blaschek, 12.6.2003
!.....................................................................

      USE alloci

      INCLUDE 'parmax.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'dat.fin'
      INCLUDE 'model.fin'       ! mit nachbar und ldir
      INCLUDE 'konv.fin'
      INCLUDE 'inv.fin'
! smatm(1:anztri,0:selmax) 0-Mitte, element über 1.,...,selmax. Kante
!.....................................................................

! PROGRAMMINTERNE PARAMETER:

      real            * 8     alfdis
      
! Hilfsvariablen 
!     Kantenlaengen und Schwerpunkte der Fl.el. nicht veraenderlich,
!     daher auch vielleicht global denkbar 
      real            * 8     dum,  
     1     edglen(selmax),      ! Laenge der Ka
     1     dist(selmax),        ! Abstand der Schwerpunkte
     1     sp(0:selmax,2)       ! Schwerpunktkoordinaten
      integer         * 4     i,j,l
      
      DOUBLE PRECISION :: ang   !angle for right smoothing alpha-estimation
      DOUBLE PRECISION :: alfxz !mixed sum for smoothing
      
!.....................................................................

      
! Rauhigkeitsmatrix auf Null setzen
      smatm(1:manz,0:selmax) = 0d0
      iflnr =0
      DO i=1,typanz
         
         IF (typ(i) <=10) THEN  ! < ? Flaechenelement
            
            DO j=1,nelanz(i)
               iflnr = iflnr + 1 
               sp(0,1)=0. ; sp(0,2)=0.
               
               DO k=1,selanz(i) ! Schwerpunkt des mittleren Elementes
                  sp(0,1)=sp(0,1)+sx(snr(nrel(iflnr,k)))
                  sp(0,2)=sp(0,2)+sy(snr(nrel(iflnr,k)))
               END DO
               
               sp(0,1)=sp(0,1)/selanz(i); sp(0,2)=sp(0,2)/selanz(i)
               DO k=1,selanz(i)
                  edglen(k)=SQRT( ( sx(snr(nrel(iflnr,k)))-  
     1                 sx(snr(nrel(iflnr,MOD(k,selanz(i))+1))) )**2+
     1                 ( sy(snr(nrel(iflnr,k)))- 
     1                 sy(snr(nrel(iflnr,MOD(k,selanz(i))+1))))**2 )
                  sp(1:selmax,1:2)=0.
c     if (iflnr==75) write(*,*) 'neighbr',nachbar(iflnr,k),ldir(iflnr,k)
                  IF (nachbar(iflnr,k)>0 .AND. ldir(iflnr,k)) THEN !existiert und soll
                     ijdum=nachbar(nachbar(iflnr,k),0)
                     
                     DO l=1,ijdum ! beruecksichtigt werden
                        sp(k,1)=sp(k,1)+
     1                       sx(snr(nrel(nachbar(iflnr,k),l)))
                        sp(k,2)=sp(k,2)+
     1                       sy(snr(nrel(nachbar(iflnr,k),l)))
                     END DO
                     
                     sp(k,1)=sp(k,1)/ijdum
                     sp(k,2)=sp(k,2)/ijdum
c     if (iflnr==75) write(*,*) 'pos',k,sp(k,1),sp(k,2)
                     dist(k)=SQRT((sp(0,1)-sp(k,1))*(sp(0,1)-sp(k,1))+ 
     1                    (sp(0,2)-sp(k,2)) * (sp(0,2)-sp(k,2)) )
                     ang=DATAN2(sp(0,2)-sp(k,2),sp(0,1)-sp(k,1)) !neu
                     alfxz=DSQRT((alfx*DCOS(ang))**2+
     1                    (alfz*DSIN(ang))**2)
!     alfxz=DSQRT((alfx)**2+(alfz)**2)
!IF (DABS(1d0-alfxz)>8d0) WRITE(*,*) 'alfxz',alfxz,DABS(1d0-alfxz),j
                     smatm(iflnr,k)=-edglen(k)/dist(k)*alfxz
c     if (iflnr==75) write(*,*) 'edganddistadf',edglen(k),dist(k)
                     smatm(iflnr,0)=smatm(iflnr,0)+
     1                    edglen(k)/dist(k)*alfxz !bis hier
!     ELSE 0-Beitrag
                  ENDIF
c     if (iflnr==75) write(*,*) ang
               END DO 
c     if (iflnr==75) write(*,*) 'Middle point ',sp(0,1),sp(0,2)
c     if (iflnr==75) write(*,*) 'elem 75: ',smatm(iflnr,:),
c     1           edglen(:),dist(:)
            END DO 
         ENDIF 
      END DO
c     write(*,*) 'Itiration ',it,', ',itr,'Smatm',(smatm(75,0))
      
      

      RETURN
      END


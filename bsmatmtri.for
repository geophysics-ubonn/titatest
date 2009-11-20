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
      
      IMPLICIT none

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
      integer         * 4     i,j,l,k,iflnr,ijdum,smaxs,ik
      REAL            * 8     edglen(selmax)   ! Kantenlaenge
      REAL            * 8     dist(selmax)     ! Abstand der Schwerpunkte
      REAL            * 8     sp(0:selmax,2)   ! Schwerpunktkoordinaten
      REAL            * 8 ang   !Winkel fuer anisotrope Glaettung
      REAL            * 8 alfxz !Anisotrope Glaettung
      
!.....................................................................
      
      smaxs=MAXVAL(selanz)

      IF (.NOT.ALLOCATED (smatm))ALLOCATE (smatm(manz,smaxs+1),STAT=errnr)
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem smatm in bsmatmtri'
         errnr = 97
         RETURN
      END IF

      smatm = 0d0               ! initialize smatm
      iflnr =0
      IF (elanz/=manz)PRINT*,'elanz/=manzSMATMTRI may be wrong'
      DO i=1,elanz

         sp(0:smaxs,:) = 0.
         
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
            

            IF (nachbar(i,k)>0) THEN !nachbar existiert 
               
               DO l=1,smaxs
                  sp(k,1) = sp(k,1) + sx(snr(nrel(nachbar(i,k),l)))
                  sp(k,2) = sp(k,2) + sy(snr(nrel(nachbar(i,k),l)))
               END DO
               
               sp(k,1) = sp(k,1)/smaxs ! schwerpunkt des nachbar elements
               sp(k,2) = sp(k,2)/smaxs
               
               dist(k) = SQRT((sp(0,1) - sp(k,1))**2 + 
     1              (sp(0,2) - sp(k,2))**2)
               
               ang = DATAN2((sp(0,2) - sp(k,2)),(sp(0,1) - sp(k,1))) !neu
               
               alfxz = DSQRT((alfx*DCOS(ang))**2 + (alfz*DSIN(ang))**2)
               
               smatm(i,k) = -edglen(k)/dist(k)*alfxz
               
               smatm(i,smaxs+1) = smatm(i,smaxs+1) + 
     1              edglen(k)/dist(k)*alfxz !Hauptdiagonale

            END IF

         END DO
      END DO

      END


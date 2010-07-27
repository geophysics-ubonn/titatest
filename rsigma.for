      subroutine rsigma(kanal,datei)

c     Unterprogramm zum Einlesen der Widerstandsverteilung aus 'datei'.

c     Andreas Kemna                                            20-Dec-1993
c     Letzte Aenderung   07-Nov-1997
      
c.....................................................................

      USE make_noise
      USE alloci,only:rnd_r,rnd_p
      USE femmod
      USE invmod

      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'sigma.fin'
      INCLUDE 'konv.fin'
      INCLUDE 'model.fin'
c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Kanalnummer
      integer         * 4     kanal

c     Datei
      character       * 80    datei

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Hilfsvariablen
      integer         * 4     i,idum,ifp1,ifp2
      real            * 8     bet,pha,eps_r,eps_p
      
c Pi
      real            * 8     pi
c.....................................................................
      pi = dacos(-1d0)
        
c     'datei' oeffnen
      fetxt = datei

      errnr = 1
      open(kanal,file=fetxt,status='old',err=999)

      errnr = 3
      IF (lnsepri) THEN
         PRINT*,''
         CALL get_unit(ifp1)
         PRINT*,'iseedpri / std. err.',iseedpri,modl_stdn
         OPEN (ifp1,FILE='tmp.mynoise_mprior_rho',STATUS='replace')
         WRITE (*,'(A)',advance='no')' Initializing noise '
         ALLOCATE (rnd_r(elanz))
         CALL Random_Init(iseedpri)
         DO i=1,elanz
            rnd_r(i) = Random_Gauss()
         END DO
         IF (.NOT.ldc) THEN
            CALL get_unit(ifp2)
            OPEN (ifp2,FILE='tmp.mynoise_mprior_phase',STATUS='replace')
            ALLOCATE (rnd_p(elanz))
            CALL Random_Init(-iseedpri)
            DO i=1,elanz
               rnd_p(i) = Random_Gauss()
            END DO
         END IF
         PRINT*,''
      END IF

c     Anzahl der Elemente (ohne Randelemente) einlesen
      read(kanal,*,end=1001,err=1000) idum

c     Ggf. Fehlermeldung
      if (idum.ne.elanz) then
         fetxt = ' '
         errnr = 47
         goto 1000
      end if

c     Betrag und Phase (in mrad) des komplexen Widerstandes einlesen
      do i=1,elanz
         read(kanal,*,end=1001,err=1000) bet,pha

         pha      = 1d-3*pha

         IF (lprior) THEN 
!     set prior model ...             
            IF (bet0 <= 0. .OR.
     1           (.not.ldc.and.dabs(pha0).gt.1d3*pi)) THEN
               fetxt = ' '
               errnr = 91
               goto 999
            END IF
            
            IF (bet > 0.) THEN  
!     TODO: meaningful phase check.. 
               IF (lnsepri) THEN
                  eps_r = 1d-2*modl_stdn * bet
                  WRITE(ifp1,'(3(G14.4,1X))',ADVANCE='no')
     1                 rnd_r(i),eps_r,bet
                  bet = bet + rnd_r(i) * eps_r
                  WRITE(ifp1,'(G14.4)')bet
                  IF (.NOT. ldc) THEN
                     eps_p = 1d-2*modl_stdn*dabs(pha)
                     WRITE(ifp2,'(3(G14.4,1X))',ADVANCE='no')
     1                    rnd_p(i),eps_p,pha
                    pha = pha + rnd_p(i) * eps_p
                    WRITE(ifp2,'(G14.4)')pha
                  END IF
               END IF
               sigma(i) = dcmplx(dcos(pha)/bet,-dsin(pha)/bet)
c    hier koennte mittelung passieren
               m0(mnr(i)) = cdlog(sigma(i))
            ELSE                
!     or let it stay at zero and take background cond
               sigma(i) = dcmplx( dcos(pha0/1d3)/bet0 ,
     1              -dsin(pha0/1d3)/bet0 )
               m0(mnr(i)) = 0.
            END IF
            
         ELSE
            
c     Ggf. Fehlermeldung
            if (bet.lt.1d-12) then
               fetxt = ' '
               errnr = 11
               goto 999
            end if
            
c     Komplexe Leitfaehigkeit bestimmen
            sigma(i) = dcmplx(dcos(pha)/bet,-dsin(pha)/bet)

         END IF
      end do

c     'datei' schliessen
      close(kanal)

      IF (lnsepri) THEN
         close (ifp1)
         IF (.NOT.ldc) close(ifp2)
      END IF

      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 999  return

 1000 close(kanal)
      return

 1001 close(kanal)
      errnr = 2
      return

      end

      subroutine rsigma(kanal,datei)

c     Unterprogramm zum Einlesen der Widerstandsverteilung aus 'datei'.

c     Andreas Kemna                                            20-Dec-1993
c     Letzte Aenderung   07-Nov-1997
      
c.....................................................................
      IMPLICIT none
      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'sigma.fin'
      INCLUDE 'konv.fin'
      INCLUDE 'model.fin'
      INCLUDE 'inv.fin'
      INCLUDE 'fem.fin'
c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Kanalnummer
      integer         * 4     kanal

c     Datei
      character       * 80    datei

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Hilfsvariablen
      integer         * 4     i,idum
      real            * 8     bet,pha
c Pi
      real            * 8     pi
c.....................................................................
      pi = dacos(-1d0)
        
c     'datei' oeffnen
      fetxt = datei

      errnr = 1
      open(kanal,file=fetxt,status='old',err=999)

      errnr = 3

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
               sigma(i) = dcmplx(dcos(pha)/bet,-dsin(pha)/bet)
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

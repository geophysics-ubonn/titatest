      subroutine bbsedc(kanal,datei)

c     Unterprogramm zur Berechnung der Summe der Sensitivitaeten aller
c     Messungen (normiert)
c     Berechnet nun coverage als Summe der Absolutbetraege..
c     ('kpotdc' als Hilfsfeld benutzt).

c     Andreas Kemna                                         02-Mar-1995
c     Letzte Aenderung                   31-Mar-2010

c.....................................................................

      USE alloci
      USE datmod
      USE invmod
      USE modelmod
      USE elemmod
      USE errmod

      IMPLICIT none


c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Kanalnummer
      integer         * 4     kanal

c     Datei
      character       * 80    datei

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Hilfsvariablen
      real            * 8     xdum,ydum,dum,dum2,dummax

c     Indexvariablen
      integer         * 4     i,j,k,l,i2

c     Aktuelle Elementnummer
      integer         * 4     iel

c     Anzahl der Knoten im aktuellen Elementtyp
      integer         * 4     nkel

c.....................................................................

c     Werte berechnen
      dum = 0d0

      do j=1,manz
         l  = j/sanz + 1
         i2 = mod(j,sanz)

         kpotdc(i2,l,1) = 0d0
         kpotdc(i2,l,2) = 0d0
         kpotdc(i2,l,3) = 0d0
         kpotdc(i2,l,4) = 0d0

         do i=1,nanz
            dum2 = SQRT(sensdc(i,j)*sensdc(i,j))*
     1           wmatd(i)*dble(wdfak(i))
            kpotdc(i2,l,1) = kpotdc(i2,l,1) + dum2
         end do

         dum = dmax1(dum,kpotdc(i2,l,1))
      end do
      dummax = dum
c     Summe der Sensitivitaeten normieren
      if (dum.gt.1d-12) then
         do j=1,manz
            l = j/sanz + 1
            i = mod(j,sanz)
            kpotdc(i,l,1) = kpotdc(i,l,1) / dum
         end do
      end if

c     Schwerpunktkoordinaten der (ggf. zusammengesetzten) Elemente bestimmen
      iel = 0

      do i=1,typanz
         if (typ(i).gt.10) goto 10

         nkel = selanz(i)

         do j=1,nelanz(i)
            iel  = iel + 1

            xdum = 0d0
            ydum = 0d0

            do k=1,nkel
               xdum = xdum + sx(snr(nrel(iel,k)))
               ydum = ydum + sy(snr(nrel(iel,k)))
            end do

            l  = mnr(iel)/sanz + 1
            i2 = mod(mnr(iel),sanz)

            kpotdc(i2,l,2) = kpotdc(i2,l,2) + xdum/dble(nkel)
            kpotdc(i2,l,3) = kpotdc(i2,l,3) + ydum/dble(nkel)
            kpotdc(i2,l,4) = kpotdc(i2,l,4) + 1d0
         end do
      end do

 10   continue

      do j=1,manz
         l = j/sanz + 1
         i = mod(j,sanz)
         kpotdc(i,l,2) = kpotdc(i,l,2) / kpotdc(i,l,4)
         kpotdc(i,l,3) = kpotdc(i,l,3) / kpotdc(i,l,4)
      end do

c     'datei' oeffnen
      fetxt = datei
      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4

c     Anzahl der Werte
      write(kanal,*,err=1000) manz

c     Koordinaten und Sensitivitaetsbetraege schreiben
c     (logarithmierter (Basis 10) normierter Betrag)
      do j=1,manz
         l   = j/sanz + 1
         i   = mod(j,sanz)
         dum = kpotdc(i,l,1)

         if (dum.gt.0d0) then
            write(kanal,*,err=1000) real(kpotdc(i,l,2)),
     1           real(kpotdc(i,l,3)),
     1           real(dlog10(dum))
         else
            write(kanal,*,err=1000) real(kpotdc(i,l,2)),
     1           real(kpotdc(i,l,3)),
     1           -1.e2
         end if
      end do

c     ak
c     Maximale Sensitivitaet schreiben
      write(kanal,*,err=1000)'Max:',dummax
c     'datei' schliessen
      close(kanal)

      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 999  return

 1000 close(kanal)
      return

      end

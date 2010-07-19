      subroutine wsens(kanal,datei)

c     Unterprogramm zum Schreiben der Sensitivitaeten aller Messungen.

c     Andreas Kemna                                            20-Apr-1995
c     Letzte Aenderung   10-Mar-2007

c.....................................................................

      USE datmod
      USE alloci

      IMPLICIT none
      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'model.fin'

c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Kanalnummer
      integer         * 4     kanal

c     Datei
      character       * 80    datei

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Hilfsfunctions
      character       * 12    intcha
      character       * 80    filpat

c     Hilfsvariablen
      integer         * 4     idum,idum2,lnanz
      character       * 80    htxt
      character       * 12    htxt2

      real            * 8     xdum,ydum,
     1     sensmax
      complex         * 16    summe
      integer         * 4     is(mmax)

c     Schwerpunktkoordinaten
      real            * 8     xs(mmax),
     1     ys(mmax)

c     Indexvariablen
      integer         * 4     i,j,k

c     Aktuelle Elementnummer
      integer         * 4     iel

c     Anzahl der Knoten im aktuellen Elementtyp
      integer         * 4     nkel

c     (Back-) Slash
      character       * 1     slash

c.....................................................................

c     Slash
      slash = '/'

c     'datei' modifizieren
      lnanz = int(log10(real(nanz)))+1
      htxt  = filpat(datei,idum2,1,slash(1:1))
      idum  = idum2+index(datei(idum2+1:80),'.')-1

      if ((idum-idum2-1).gt.(8-lnanz)) then
         fetxt = datei
         errnr = 15
         goto 999
      end if

c     Schwerpunktkoordinaten der (ggf. zusammengesetzten) Elemente bestimmen
      do i=1,manz
         xs(i) = 0d0
         ys(i) = 0d0
         is(i) = 0
      end do

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

            xs(mnr(iel)) = xs(mnr(iel)) + xdum / dble(nkel)
            ys(mnr(iel)) = ys(mnr(iel)) + ydum / dble(nkel)

            is(mnr(iel)) = is(mnr(iel)) + 1
         end do
      end do

 10   continue

      do i=1,manz
         xs(i) = xs(i) / dble(is(i))
         ys(i) = ys(i) / dble(is(i))
      end do

c     Ausgabeschleife
      do i=1,nanz

c     'datei' bestimmen
         htxt2 = intcha(i,lnanz)
         htxt  = datei(1:idum)//htxt2(1:lnanz)//datei(idum+1:idum+4)

c     'datei' oeffnen
         fetxt = htxt
         errnr = 1
         open(kanal,file=fetxt,status='replace',err=999)
         errnr = 4

c     Summe der Sensitivitaeten berechnen und ausgeben
c     (Betrag und Phase (in mrad))
c     ak            sensmax = 0d0
         summe   = dcmplx(0d0)

         do j=1,manz
c     ak                sensmax = dmax1(sensmax,cdabs(sens(i,j)))
            summe   = summe + sens(i,j)
         end do

         if (cdabs(summe).gt.0d0) then
            write(kanal,*,err=1000) real(cdabs(summe)),
     1           real(1d3*datan2(dimag(summe),dble(summe)))
         else
            write(kanal,*,err=1000) 0.,0.
         end if

c     Koordinaten und Sensitivitaeten schreiben
c     (Real- und Imaginaerteil)
         do j=1,manz
            write(kanal,*,err=1000)
     1           real(xs(j)),real(ys(j)),
     1           real(dble(sens(i,j))),real(dimag(sens(i,j)))
         end do

c     'datei' schliessen
         close(kanal)
      end do

      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 999  return

 1000 close(kanal)
      return

      end

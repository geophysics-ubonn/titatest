      subroutine wsens(kanal,datei)

!!!$     Unterprogramm zum Schreiben der Sensitivitaeten aller Messungen.

!!!$     Andreas Kemna                                            20-Apr-1995
!!!$     Letzte Aenderung   10-Mar-2007

!!!$.....................................................................

      USE datmod
      USE alloci
      USE modelmod
      USE elemmod
      USE errmod

      IMPLICIT none
      

!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Kanalnummer
      integer         * 4     kanal

!!!$     Datei
      character       * 80    datei

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsfunctions
      character       * 12    intcha
      character       * 80    filpat

!!!$     Hilfsvariablen
      integer         * 4     idum,idum2,lnanz
      character       * 80    htxt
      character       * 12    htxt2

      real            * 8     xdum,ydum,
     1     sensmax
      complex         * 16    summe
!!!$     Indexvariablen
      integer         * 4     i,j,k

!!!$     Aktuelle Elementnummer
      integer         * 4     iel

!!!$     Anzahl der Knoten im aktuellen Elementtyp
      integer         * 4     nkel

!!!$     (Back-) Slash
      character       * 1     slash

!!!$.....................................................................

!!!$     Slash
      slash = '/'

!!!$     'datei' modifizieren
      lnanz = int(log10(real(nanz)))+1
      htxt  = filpat(datei,idum2,1,slash(1:1))
      idum  = idum2+index(datei(idum2+1:80),'.')-1

      if ((idum-idum2-1).gt.(8-lnanz)) then
         fetxt = datei
         errnr = 15
         goto 999
      end if

!!!$     Ausgabeschleife
      do i=1,nanz

!!!$     'datei' bestimmen
         htxt2 = intcha(i,lnanz)
         htxt  = datei(1:idum)//htxt2(1:lnanz)//datei(idum+1:idum+4)

!!!$     'datei' oeffnen
         fetxt = htxt
         errnr = 1
         open(kanal,file=fetxt,status='replace',err=999)
         errnr = 4

!!!$     Summe der Sensitivitaeten berechnen und ausgeben
!!!$     (Betrag und Phase (in mrad))
!!!$     ak            sensmax = 0d0
         summe   = dcmplx(0d0)

         do j=1,manz
!!!$     ak                sensmax = dmax1(sensmax,cdabs(sens(i,j)))
            summe   = summe + sens(i,j)
         end do

         if (cdabs(summe).gt.0d0) then
            write(kanal,*,err=1000) real(cdabs(summe)),
     1           real(1d3*datan2(dimag(summe),dble(summe)))
         else
            write(kanal,*,err=1000) 0.,0.
         end if

!!!$     Koordinaten und Sensitivitaeten schreiben
!!!$     (Real- und Imaginaerteil)
         do j=1,manz
            write(kanal,*,err=1000)
     1           real(espx(j)),real(espy(j)),
     1           real(dble(sens(i,j))),real(dimag(sens(i,j)))
         end do

!!!$     'datei' schliessen
         close(kanal)
      end do

      errnr = 0
      return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

 999  return

 1000 close(kanal)
      return

      end

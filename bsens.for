      subroutine bsens()

c     Unterprogramm zur Berechnung der Sensitivitaeten.

c     Andreas Kemna                                            09-Apr-1995
c     Letzte Aenderung   04-Dez-1998

c.....................................................................

      USE alloci
      USE femmod
      USE datmod

      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'electr.fin'
      INCLUDE 'waven.fin'
      INCLUDE 'model.fin'

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Aktuelle Elementnummer
      integer         * 4     iel

c     Aktueller Elementtyp
      integer         * 4     ntyp

c     Anzahl der Knoten im aktuellen Elementtyp
      integer         * 4     nkel

c     Elektrodennummern
      integer         * 4     elec1,elec2,
     1     elec3,elec4

c     Beitraege zur Superposition
      complex         * 16    sup(4)

c     Indexvariablen
      integer         * 4     ityp,jnel,mi,mj,
     1     imn,imax,imin
      integer         * 4     i,j,k

c     Hilfsfeld
      complex         * 16    hsens(kwnmax)

c     Hilfsvariablen
      integer         * 4     nzp,nnp
      complex         * 16    dum

c     Pi
      real            * 8     pi

c.....................................................................

      pi = dacos(-1d0)

c     Sensitivitaetenfeld auf Null setzen
      do i=1,nanz
         do j=1,manz
            sens(i,j) = dcmplx(0d0)
         end do
      end do

c     Messwert hochzaehlen
      do i=1,nanz
         iel = 0

c     Kontrollausgabe
         write(*,'(''+'',a,i4)',ADVANCE='no')ACHAR(13)//
     1        ' Calculating Sensitivities : Measurement ',i

c     Stromelektroden bestimmen
         elec1 = mod(strnr(i),10000)
         elec2 = (strnr(i)-elec1)/10000

c     Messelektroden bestimmen
         elec3 = mod(vnr(i),10000)
         elec4 = (vnr(i)-elec3)/10000

c     Beitraege zur Superposition auf Null setzen
         do j=1,4
            sup(j) = dcmplx(0d0)
         end do

         do ityp=1,typanz
            ntyp = typ(ityp)
            nkel = selanz(ityp)

c     Ggf. zu neuem Messwert springen
            if (ntyp.gt.10) goto 10

            do jnel=1,nelanz(ityp)

c     Elementnummer hochzaehlen
               iel = iel + 1

c     SENSITIVITAETEN BERECHNEN
               do k=1,kwnanz
                  hsens(k) = dcmplx(0d0)

c     Knoten des aktuellen Elements hochzaehlen
                  do mi=1,nkel
                     nzp = nrel(iel,mi)

                     do mj=1,nkel
                        nnp  = nrel(iel,mj)
                        imax = max0(mi,mj)
                        imin = min0(mi,mj)
                        imn  = imax*(imax-1)/2+imin

c     Beitraege nach "Reziprozitaetsmethode" gewichtet aufaddieren und
c     superponieren
c     (beachte: 'volt = pot(elec4) - pot(elec3)' ,
c     '+I' bei 'elec2', '-I' bei 'elec1' )
                        if (elec1.gt.0) sup(1) = kpot(nnp,elec1,k)
                        if (elec2.gt.0) sup(2) = kpot(nnp,elec2,k)
                        if (elec3.gt.0) sup(3) = kpot(nzp,elec3,k)
                        if (elec4.gt.0) sup(4) = kpot(nzp,elec4,k)

c     ACHTUNG: Bei grossen Quellabstaenden UNDERFLOW moeglich, da einzelnen
c     Potentiale sehr klein (vor allem bei grossen Wellenzahlen)!
c     -> mittels Compiler-Einstellung auf Null setzen!
c     MsDev5.0: "/fpe:3 /check:underflow" -> "/fpe:0"
                        dum      = (sup(2)-sup(1)) * (sup(4)-sup(3))
                        hsens(k) = hsens(k) + dcmplx(elbg(iel,imn,k))
     1                       * dum
                     end do
                  end do
               end do

c     GGF. RUECKTRANSFORMATION
               if (swrtr.eq.0) then

                  dum = hsens(1) * dcmplx(strom(i))

               else

                  dum = dcmplx(0d0)

                  do k=1,kwnanz
                     dum = dum + hsens(k)*dcmplx(kwnwi(k))
                  end do

                  dum = dum * dcmplx(strom(i)/pi)

               end if

c     Sensitivitaeten aufaddieren
               sens(i,mnr(iel)) = sens(i,mnr(iel)) - dum

            end do
         end do

 10      continue
      end do

      return
      end

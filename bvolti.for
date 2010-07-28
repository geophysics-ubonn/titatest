      subroutine bvolti()

c     Unterprogramm zur Berechnung der modellierten Spannungswerte
c     (beachte: Potentialwerte wurden fuer Einheitsstrom berechnet).

c     Andreas Kemna                                            03-Sep-1994
c     Letzte Aenderung   13-Nov-1997
      
c.....................................................................

      USE alloci
      USE femmod
      USE datmod
      USE electrmod

      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Indexvariable
      integer         * 4     i

c     Elektrodennummern
      integer         * 4     elec1,elec2,
     1     elec3,elec4

c     Hilfsvariablen
      complex         * 16    dum1,dum2,dum3,dum4
      real            * 8     dum1dc,dum2dc,dum3dc,dum4dc

c.....................................................................

      do i=1,nanz

c     Stromelektroden bestimmen
         elec1 = mod(strnr(i),10000)
         elec2 = (strnr(i)-elec1)/10000

c     Messelektroden bestimmen
         elec3 = mod(vnr(i),10000)
         elec4 = (vnr(i)-elec3)/10000

c     Spannungswert berechnen (Superposition)
c     (beachte: Faktoren '.../2d0' (-> Potentialwerte fuer Einheitsstrom)
c     und '...*2d0' (-> Ruecktransformation) kuerzen sich weg !)
         if (ldc) then
            dum1dc = dble(min0(elec4,1)*min0(elec2,1))
     1           *hpotdc(enr(max0(elec4,1)),max0(elec2,1))
            dum2dc = dble(min0(elec4,1)*min0(elec1,1))
     1           *hpotdc(enr(max0(elec4,1)),max0(elec1,1))
            dum3dc = dble(min0(elec3,1)*min0(elec2,1))
     1           *hpotdc(enr(max0(elec3,1)),max0(elec2,1))
            dum4dc = dble(min0(elec3,1)*min0(elec1,1))
     1           *hpotdc(enr(max0(elec3,1)),max0(elec1,1))

            volt(i) = dcmplx((dum1dc-dum2dc) - (dum3dc-dum4dc))
         else
            dum1 = dcmplx(min0(elec4,1)*min0(elec2,1))
     1           *hpot(enr(max0(elec4,1)),max0(elec2,1))
            dum2 = dcmplx(min0(elec4,1)*min0(elec1,1))
     1           *hpot(enr(max0(elec4,1)),max0(elec1,1))
            dum3 = dcmplx(min0(elec3,1)*min0(elec2,1))
     1           *hpot(enr(max0(elec3,1)),max0(elec2,1))
            dum4 = dcmplx(min0(elec3,1)*min0(elec1,1))
     1           *hpot(enr(max0(elec3,1)),max0(elec1,1))

            volt(i) = (dum1-dum2) - (dum3-dum4)
         end if

         if (cdabs(volt(i)).eq.0d0) then
            fetxt = 'index '
            write(fetxt(7:10),'(i4)') i
            errnr = 82
            goto 1000
         end if

c     Werte logarithmieren
         sigmaa(i) = -cdlog(volt(i))
      end do

      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 1000 return

      end

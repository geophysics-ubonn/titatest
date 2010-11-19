      subroutine bvolt()

c     Unterprogramm zur Berechnung der Spannungs- und scheinbaren
c     Widerstandswerte (beachte: Potentialwerte wurden fuer Einheitsstrom
c     berechnet).

c     Andreas Kemna                                            03-Sep-1994
c     Letzte Aenderung   06-Nov-1997

c.....................................................................

      USE alloci
      USE datmod
      USE electrmod
      USE errmod

      IMPLICIT none


c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Indexvariable
      INTEGER (KIND = 4)  ::     i,j

c     Elektrodennummern
      INTEGER (KIND = 4)  ::     elec1,elec2,
     1     elec3,elec4

c     Hilfsvariablen
      COMPLEX (KIND(0D0)) ::    dum1,dum2,dum3,dum4

c.....................................................................
      j=0
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
         dum1 = dcmplx(min0(elec4,1)*min0(elec2,1))
     1        *hpot(enr(max0(elec4,1)),max0(elec2,1))
         dum2 = dcmplx(min0(elec4,1)*min0(elec1,1))
     1        *hpot(enr(max0(elec4,1)),max0(elec1,1))
         dum3 = dcmplx(min0(elec3,1)*min0(elec2,1))
     1        *hpot(enr(max0(elec3,1)),max0(elec2,1))
         dum4 = dcmplx(min0(elec3,1)*min0(elec1,1))
     1        *hpot(enr(max0(elec3,1)),max0(elec1,1))

         volt(i) = (dum1-dum2) - (dum3-dum4)

         if (cdabs(volt(i)).eq.0d0) then
            j=j+1
c     ak
            write(*,'(A,I8,A,I8)',ADVANCE='no')ACHAR(13)//
     1           ' --- Messpannung',i,' ist null ',j
c     RM               fetxt = ' '
c     RM               errnr = 82
c     RM               goto 1000
         end if

c     Scheinbaren Widerstandswert berechnen
c     ak            sigmaa(i) = dcmplx(kfak(i)) * volt(i)
         sigmaa(i) = volt(i)
         sigmaa(i) = dcmplx(kfak(i)) * volt(i)
      end do
      IF (j/=0) WRITE (*,'(/A,I8,A)')
     1     ' Vorsicht es wurde',j,'mal keine Messpannung gemessen'
      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 1000 return

      end

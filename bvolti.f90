subroutine bvolti()

!!!$     Unterprogramm zur Berechnung der modellierten Spannungswerte
!!!$     (beachte: Potentialwerte wurden fuer Einheitsstrom berechnet).

!!!$     Andreas Kemna                                            03-Sep-1994
!!!$     Letzte Aenderung   13-Nov-1997

!!!$.....................................................................

  USE alloci
  USE femmod
  USE datmod
  USE electrmod
  USE errmod

  IMPLICIT none


!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Indexvariable
  INTEGER (KIND = 4)  ::     i

!!!$     Elektrodennummern
  INTEGER (KIND = 4)  ::     elec1,elec2,elec3,elec4

!!!$     Hilfsvariablen
  COMPLEX (KIND(0D0)) ::    dum1,dum2,dum3,dum4
  REAL (KIND(0D0))    ::     dum1dc,dum2dc,dum3dc,dum4dc

!!!$.....................................................................

  do i=1,nanz

!!!$     Stromelektroden bestimmen
     elec1 = mod(strnr(i),10000)
     elec2 = (strnr(i)-elec1)/10000

!!!$     Messelektroden bestimmen
     elec3 = mod(vnr(i),10000)
     elec4 = (vnr(i)-elec3)/10000

!!!$     Spannungswert berechnen (Superposition)
!!!$     (beachte: Faktoren '.../2d0' (-> Potentialwerte fuer Einheitsstrom)
!!!$     und '...*2d0' (-> Ruecktransformation) kuerzen sich weg !)
     if (ldc) then
        dum1dc = dble(min0(elec4,1)*min0(elec2,1)) &
             *hpotdc(enr(max0(elec4,1)),max0(elec2,1))
        dum2dc = dble(min0(elec4,1)*min0(elec1,1)) &
             *hpotdc(enr(max0(elec4,1)),max0(elec1,1))
        dum3dc = dble(min0(elec3,1)*min0(elec2,1)) &
             *hpotdc(enr(max0(elec3,1)),max0(elec2,1))
        dum4dc = dble(min0(elec3,1)*min0(elec1,1)) &
             *hpotdc(enr(max0(elec3,1)),max0(elec1,1))

        volt(i) = dcmplx((dum1dc-dum2dc) - (dum3dc-dum4dc))
     else
        dum1 = dcmplx(min0(elec4,1)*min0(elec2,1)) &
             *hpot(enr(max0(elec4,1)),max0(elec2,1))
        dum2 = dcmplx(min0(elec4,1)*min0(elec1,1)) &
             *hpot(enr(max0(elec4,1)),max0(elec1,1))
        dum3 = dcmplx(min0(elec3,1)*min0(elec2,1)) &
             *hpot(enr(max0(elec3,1)),max0(elec2,1))
        dum4 = dcmplx(min0(elec3,1)*min0(elec1,1)) &
             *hpot(enr(max0(elec3,1)),max0(elec1,1))

        volt(i) = (dum1-dum2) - (dum3-dum4)
     end if

     if (cdabs(volt(i)).eq.0d0) then
        write(fetxt,*)'index',i
        print*,hpotdc(enr(i),:),dum1dc,dum2dc,dum3dc,dum4dc
        errnr = 82
        goto 1000
     end if

!!!$     Werte logarithmieren
     sigmaa(i) = -cdlog(volt(i))
  end do

  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

1000 return

end subroutine bvolti

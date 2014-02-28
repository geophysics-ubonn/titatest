SUBROUTINE bvolti()

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

  IMPLICIT NONE


!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Indexvariable
  INTEGER (KIND = 4)  ::     i

!!!$     Elektrodennummern
  INTEGER (KIND = 4)  ::     elec1,elec2,elec3,elec4

!!!$     Hilfsvariablen
  COMPLEX (prec) ::    dum1,dum2,dum3,dum4
  REAL (prec)    ::     dum1dc,dum2dc,dum3dc,dum4dc

!!!$.....................................................................

  DO i=1,nanz

!!!$     Stromelektroden bestimmen
     elec1 = MOD(strnr(i),10000)
     elec2 = (strnr(i)-elec1)/10000

!!!$     Messelektroden bestimmen
     elec3 = MOD(vnr(i),10000)
     elec4 = (vnr(i)-elec3)/10000

!!!$     Spannungswert berechnen (Superposition)
!!!$     (beachte: Faktoren '.../2d0' (-> Potentialwerte fuer Einheitsstrom)
!!!$     und '...*2d0' (-> Ruecktransformation) kuerzen sich weg !)
     IF (ldc) THEN
        dum1dc = real(min0(elec4,1)*min0(elec2,1)) &
             *hpotdc(enr(max0(elec4,1)),max0(elec2,1))
        dum2dc = real(min0(elec4,1)*min0(elec1,1)) &
             *hpotdc(enr(max0(elec4,1)),max0(elec1,1))
        dum3dc = real(min0(elec3,1)*min0(elec2,1)) &
             *hpotdc(enr(max0(elec3,1)),max0(elec2,1))
        dum4dc = real(min0(elec3,1)*min0(elec1,1)) &
             *hpotdc(enr(max0(elec3,1)),max0(elec1,1))

        volt(i) = CMPLX((dum1dc-dum2dc) - (dum3dc-dum4dc))
     ELSE
        dum1 = CMPLX(min0(elec4,1)*min0(elec2,1)) &
             *hpot(enr(max0(elec4,1)),max0(elec2,1))
        dum2 = CMPLX(min0(elec4,1)*min0(elec1,1)) &
             *hpot(enr(max0(elec4,1)),max0(elec1,1))
        dum3 = CMPLX(min0(elec3,1)*min0(elec2,1)) &
             *hpot(enr(max0(elec3,1)),max0(elec2,1))
        dum4 = CMPLX(min0(elec3,1)*min0(elec1,1)) &
             *hpot(enr(max0(elec3,1)),max0(elec1,1))

        volt(i) = (dum1-dum2) - (dum3-dum4)
     END IF

     IF (ABS(volt(i)).EQ.0d0) THEN
        WRITE(fetxt,*)'index',i
        PRINT*,'bvolti: measured voltage = 0 for config',i
        errnr = 82
        GOTO 1000
     END IF

!!!$     Werte logarithmieren
     sigmaa(i) = -LOG(volt(i))
  END DO

  errnr = 0
  RETURN

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

1000 RETURN

END SUBROUTINE bvolti

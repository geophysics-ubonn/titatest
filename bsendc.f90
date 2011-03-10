subroutine bsendc()

!!!$     Unterprogramm zur Berechnung der Sensitivitaeten.

!!!$     Andreas Kemna                                            09-Apr-1995
!!!$     Letzte Aenderung   07-Mar-2003

!!!$.....................................................................

  USE alloci
  USE femmod
  USE datmod
  USE sigmamod
  USE electrmod
  USE modelmod
  USE elemmod
  USE wavenmod
  USE errmod
  USE konvmod , ONLY : lverb
  IMPLICIT none

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Aktuelle Elementnummer
  INTEGER (KIND = 4)  :: iel

!!!$     Aktueller Elementtyp
  INTEGER (KIND = 4)  :: ntyp

!!!$     Anzahl der Knoten im aktuellen Elementtyp
  INTEGER (KIND = 4)  :: nkel

!!!$     Elektrodennummern
  INTEGER (KIND = 4)  :: elec1,elec2,elec3,elec4

!!!$     Beitraege zur Superposition
  REAL (KIND(0D0))    ::    sup(4)

!!!$     Indexvariablen
  INTEGER (KIND = 4)  :: ityp,jnel,mi,mj,imn,imax,imin
  INTEGER (KIND = 4)  :: i,j,k,count

!!!$     Hilfsfeld
  REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: hsens

!!!$     Hilfsvariablen
  INTEGER (KIND = 4)  :: nzp,nnp
  REAL (KIND(0D0))    :: dum

!!!$     Pi
  REAL (KIND(0D0))    :: pi

!!!$.....................................................................

  pi = dacos(-1d0)

  !     get memory for hsens
  ALLOCATE (hsens(kwnanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation hsens'
     errnr = 94
     RETURN
  END IF

!!!$     Sensitivitaetenfeld auf Null setzen
  sensdc = 0D0
  count  = 0
  !$OMP PARALLEL DEFAULT (SHARED) &
  !$OMP FIRSTPRIVATE (hsens) &
  !$OMP PRIVATE(iel,elec1,elec2,elec3,elec4,sup,ntyp,jnel,nkel,nzp,nnp,imax,dum)
  !$OMP DO
!!!$     Messwert hochzaehlen
  do i=1,nanz

     !$OMP ATOMIC
     count = count + 1

     IF (lverb) write(*,'(a,t70,F6.2,A)',advance='no')ACHAR(13)//&
          'Sensitivity/ ',REAL( count * (100./nanz)),'%'
     iel = 0

!!!$     Stromelektroden bestimmen
     elec1 = mod(strnr(i),10000)
     elec2 = (strnr(i)-elec1)/10000

!!!$     Messelektroden bestimmen
     elec3 = mod(vnr(i),10000)
     elec4 = (vnr(i)-elec3)/10000

!!!$     Beitraege zur Superposition auf Null setzen
     do j=1,4
        sup(j) = 0d0
     end do

     do ityp=1,typanz
        ntyp = typ(ityp)
        nkel = selanz(ityp)

!!!$     Ggf. zu neuem Messwert springen
        if (ntyp.gt.10) CYCLE

        do jnel=1,nelanz(ityp)

!!!$     Elementnummer hochzaehlen
           iel = iel + 1

!!!$     SENSITIVITAETEN BERECHNEN
           do k=1,kwnanz
              hsens(k) = 0d0

!!!$     Knoten des aktuellen Elements hochzaehlen
              do mi=1,nkel
                 nzp = nrel(iel,mi)

                 do mj=1,nkel
                    nnp  = nrel(iel,mj)
                    imax = max0(mi,mj)
                    imin = min0(mi,mj)
                    imn  = imax*(imax-1)/2+imin

!!!$     Beitraege nach "Reziprozitaetsmethode" gewichtet aufaddieren und
!!!$     superponieren
!!!$     (beachte: 'volt = pot(elec4) - pot(elec3)' ,
!!!$     '+I' bei 'elec2', '-I' bei 'elec1' )
                    if (elec1.gt.0) sup(1) = kpotdc(nnp,elec1,k)
                    if (elec2.gt.0) sup(2) = kpotdc(nnp,elec2,k)
                    if (elec3.gt.0) sup(3) = kpotdc(nzp,elec3,k)
                    if (elec4.gt.0) sup(4) = kpotdc(nzp,elec4,k)

!!!$     ACHTUNG: Bei grossen Quellabstaenden UNDERFLOW moeglich, da einzelnen
!!!$     Potentiale sehr klein (vor allem bei grossen Wellenzahlen)!
!!!$     -> mittels Compiler-Einstellung auf Null setzen!
!!!$     MsDev5.0: "/fpe:3 /check:underflow" -> "/fpe:0"
                    dum      = (sup(2)-sup(1)) * (sup(4)-sup(3))
                    hsens(k) = hsens(k) + elbg(iel,imn,k) * dum
                 end do
              end do
           end do

!!!$     GGF. RUECKTRANSFORMATION
           if (swrtr.eq.0) then

              dum = hsens(1)

           else

              dum = 0d0

              do k=1,kwnanz
                 dum = dum + hsens(k)*kwnwi(k)
              end do

              dum = dum / pi

           end if

           sensdc(i,mnr(iel)) = sensdc(i,mnr(iel)) + dum * &
                dble(sigma(iel)/volt(i))

!!!$     ak BAW-Tank
!!!$     ak                if (mnr(iel).le.14*58) sensdc(i,mnr(iel))=0d0

        end do ! jnel=1,nelanz(i)
     end do ! ityp=1,typanz
  end do ! i=1,nanz
  !$OMP END PARALLEL

  DEALLOCATE (hsens)
  
end subroutine bsendc

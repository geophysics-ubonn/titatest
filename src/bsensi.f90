SUBROUTINE bsensi(tictoc)

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
  USE konvmod,ONLY: lverb
  USE ompmod
  USE tic_toc

  IMPLICIT NONE

!!!$.....................................................................
  LOGICAL :: tictoc
!!!$     PROGRAMMINTERNE PARAMETER:
  INTEGER,PARAMETER :: ntd=2
!!!$     Aktuelle Elementnummer
  INTEGER (KIND = 4)  ::     iel

!!!$     Aktueller Elementtyp
  INTEGER (KIND = 4)  ::     ntyp

!!!$     Anzahl der Knoten im aktuellen Elementtyp
  INTEGER (KIND = 4)  ::     nkel

!!!$     Elektrodennummern
  INTEGER (KIND = 4)  ::     elec1,elec2,elec3,elec4

!!!$     Beitraege zur Superposition
  COMPLEX (prec) ::    sup(4)

!!!$     Indexvariablen
  INTEGER (KIND = 4)  ::     ityp,jnel,mi,mj,imn,imax,imin
  INTEGER (KIND = 4)  ::     i,j,k,icount

!!!$     Hilfsfeld
  COMPLEX(prec),DIMENSION(:),ALLOCATABLE :: hsens

!!!$     Hilfsvariablen
  INTEGER (KIND = 4)  ::     nzp,nnp,c1
  COMPLEX (prec) ::    dum

!!!$     Pi
  REAL (prec)    ::     pi

!!!$.....................................................................

  pi = dacos(-1d0)
  !     get memory for hsens
  ALLOCATE (hsens(kwnanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation hsens'
     errnr = 94
     RETURN
  END IF

  IF (tictoc) CALL tic(c1)
!!!$     Sensitivitaetenfeld auf Null setzen
  sens = dcmplx(0d0)

  icount = 0
  !$OMP PARALLEL DEFAULT (none) &
  !$OMP FIRSTPRIVATE (hsens) &
  !$OMP SHARED (icount,strnr,vnr,typanz,typ,selanz,kwnanz,lverb,sigma,&
  !$OMP nrel,kpot,elbg,strom,kwnwi,pi,mnr,nanz,swrtr,nelanz,volt,sens) &
  !$OMP PRIVATE(iel,elec1,elec2,elec3,elec4,sup,imin,imn,&
  !$OMP ntyp,jnel,nkel,nzp,nnp,imax,dum)
  !$OMP DO SCHEDULE (GUIDED,CHUNK_0)
!!!$     Messwert hochzaehlen
  DO i=1,nanz
     iel = 0

     !$OMP ATOMIC
     icount = icount + 1

     IF (lverb) WRITE(*,'(a,t50,F6.2,A)',advance='no')ACHAR(13)//&
          'Sensitivity/ ',REAL(icount)/REAL(nanz) * 100.,'%'

!!!$     Stromelektroden bestimmen
     elec1 = MOD(strnr(i),10000)
     elec2 = (strnr(i)-elec1)/10000

!!!$     Messelektroden bestimmen
     elec3 = MOD(vnr(i),10000)
     elec4 = (vnr(i)-elec3)/10000

!!!$     Beitraege zur Superposition auf Null setzen
     DO j=1,4
        sup(j) = dcmplx(0d0)
     END DO

     DO ityp=1,typanz
        ntyp = typ(ityp)
        nkel = selanz(ityp)

!!!$     Ggf. zu neuem Messwert springen
        IF (ntyp.GT.10) CYCLE

        DO jnel=1,nelanz(ityp)

!!!$     Elementnummer hochzaehlen
           iel = iel + 1

!!!$     SENSITIVITAETEN BERECHNEN
           DO k=1,kwnanz
              hsens(k) = dcmplx(0d0)

!!!$     Knoten des aktuellen Elements hochzaehlen
imn = 0
              DO mi=1,nkel
                 nzp = nrel(iel,mi)

                 DO mj=1,nkel
                    nnp  = nrel(iel,mj)
!                    imax = max0(mi,mj)
!                    imin = min0(mi,mj)
!                    imn  = imax*(imax-1)/2+imin
imn = imn+1
!!!$     Beitraege nach "Reziprozitaetsmethode" gewichtet aufaddieren und
!!!$     superponieren
!!!$     (beachte: 'volt = pot(elec4) - pot(elec3)' ,
!!!$     '+I' bei 'elec2', '-I' bei 'elec1' )
                    IF (elec1.GT.0) sup(1) = kpot(nnp,elec1,k)
                    IF (elec2.GT.0) sup(2) = kpot(nnp,elec2,k)
                    IF (elec3.GT.0) sup(3) = kpot(nzp,elec3,k)
                    IF (elec4.GT.0) sup(4) = kpot(nzp,elec4,k)

!!!$     ACHTUNG: Bei grossen Quellabstaenden UNDERFLOW moeglich, da einzelnen
!!!$     Potentiale sehr klein (vor allem bei grossen Wellenzahlen)!
!!!$     -> mittels Compiler-Einstellung auf Null setzen!
!!!$     MsDev5.0: "/fpe:3 /check:underflow" -> "/fpe:0"
                    dum      = (sup(2)-sup(1)) * (sup(4)-sup(3))
                    hsens(k) = hsens(k) + dcmplx(elbg(iel,imn,k))* dum

                 END DO
              END DO
           END DO
!!!$     GGF. RUECKTRANSFORMATION
           IF (swrtr.EQ.0) THEN

              dum = hsens(1)

           ELSE

              dum = dcmplx(0d0)

              DO k=1,kwnanz
                 dum = dum + hsens(k)*dcmplx(kwnwi(k))
              END DO

              dum = dum / dcmplx(pi)

           END IF
!!!$     hier koennte auch eine mittelung passieren
           sens(i,mnr(iel)) = sens(i,mnr(iel)) + dum * sigma(iel)/volt(i)
        END DO
     END DO
  END DO
  !$OMP END PARALLEL
  fetxt = 'bsensi::'
  IF (tictoc) CALL toc(c1,fetxt)
  DEALLOCATE (hsens)

END SUBROUTINE bsensi

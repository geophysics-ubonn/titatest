MODULE bsmatm_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!$ Collection of subroutines to set the Regularization matrix 
!!!$ (smatm) on different purposes.
!!!$ 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!$ Copyright by Andreas Kemna 2010
!!!$
!!!$ Created by Roland Martin               30-Jul-2010
!!!$
!!!$ Last changed       RM                  Jul-2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE tic_toc ! counts calculation time
  USE alloci , ONLY : sens,sensdc,smatm,nachbar
  USE femmod , ONLY : fak,ldc
  USE elemmod, ONLY : smaxs,sx,sy,espx,espy,nrel,snr,elanz
  USE invmod , ONLY : lip,par,wmatd,wdfak
  USE errmod , ONLY : errnr,fetxt
  USE konvmod , ONLY : ltri,lgauss,lam,nx,nz,alfx,alfz,betamgs
  USE modelmod , ONLY : manz
  USE datmod , ONLY : nanz
  USE errmod, ONLY : errnr,fetxt
  USE sigmamod , ONLY : sigma
  USE variomodel 
  USE pathmod

  IMPLICIT none

  REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE,PRIVATE :: csens 

  PUBLIC :: bsmatm 
!!!$ controls which regularization is to apply


  PRIVATE :: bcsens
!!$ calculates daigonal of A^TC_d^{-1}A

  PRIVATE :: bsmatmreg
!!!$ sets smatm for smooth regularization on regular grids  
  PRIVATE :: bsmatmtri
!!!$ .. same but for unstructured grids (recommended)
  PRIVATE :: bsmatmlma
!!!$ Levenberg and Levenber-Marquardt Regularization
  PRIVATE :: bsmatmmgs
!!!$ Minimum gradient support regu
  PRIVATE :: bsmatmtv
!!!$ Total variance regu (alpha)
  PRIVATE :: bsmatmsto
!!!$ Stochastical regularization

CONTAINS  

  SUBROUTINE bsmatm(it)
!!!$
!!!$ This sub is the control unit of the smatm calculation
!!!$
    INTEGER (KIND = 4 ),INTENT(IN) :: it
    INTEGER (KIND = 4 )   :: c1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!$    get time
<<<<<<< HEAD
    errnr = 2
=======
>>>>>>> bdb481e37f3eee805916831adfae93cd60a9affc

    CALL TIC(c1)

    IF (.NOT.ALLOCATED(csens)) ALLOCATE (csens(manz),STAT=errnr)
    IF (errnr/=0) THEN
       fetxt = 'Allocation problem csens in bsmatmmgs'
       WRITE (*,'(/a/)')TRIM(fetxt)
       errnr = 97
       RETURN
    END IF
    
    IF (it == 1) WRITE (*,'(/a)',ADVANCE='no')'Regularization::'

    IF (ltri == 0) THEN
       
       WRITE (*,'(a)')' Rectangular smooth'
       CALL bsmatmreg
       
    ELSE IF (ltri == 1.OR.ltri == 2) THEN
       
       WRITE (*,'(a)')' Triangular smooth'
       CALL bsmatmtri
       
    ELSE IF (ltri == 3.OR.ltri==4) THEN
       
       IF (it == 1) THEN
          IF (ltri == 3) WRITE (*,'(a)') &
               'Levenberg damping (LA)'
          IF (ltri == 4) WRITE (*,'(a)') &
               'Marquardt-Levenberg damping (LMA)'
       ELSE
          WRITE (*,'(a)')' Updating damping'
       END IF
       
       CALL bsmatmlma
       
    ELSE IF (ltri > 4 .AND. ltri < 10) THEN
       
       IF (it == 1) THEN
          IF (ltri == 5) WRITE (*,'(a)') &
               ' Triangular pure MGS (beta)'
          IF (ltri == 6.OR.ltri == 8) WRITE (*,'(a)') &
               ' Triangular sens weighted MGS (beta)'
          IF (ltri == 7.OR.ltri == 9) WRITE (*,'(a)') &
               ' Triangular sens weighted MGS mean (beta)'
       ELSE
          WRITE (*,'(a)')' Updating MGS functional'
       END IF
       
       CALL bsmatmmgs
       
    ELSE IF (ltri == 10) THEN
       
       WRITE (*,'(a)')' Triangular Total variance (alpha)'
       CALL bsmatmtv

       
    ELSE IF (ltri == 15) THEN
       
       WRITE (*,'(a)')' Triangular Stochastic (beta)'
       CALL bsmatmsto

       IF (errnr /= 0) STOP

    ELSE
       
       WRITE (*,'(a)')' Error:: '// &
            'Regularization can just be '//&
            '0,1,3,4,5,6,7,8,9,10 or 15'
       RETURN
       
    END IF
    
    IF (ALLOCATED(csens)) DEALLOCATE (csens)

    fetxt = 'C_m^ calculation time'
    CALL TOC(c1,fetxt)

  END SUBROUTINE bsmatm

!!! 
  SUBROUTINE bcsens (csensmax,csensavg)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!$ This subroutine calculates the squared coverage or diag{A^TC_d^{-1}A}!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER :: i,j
    REAL(KIND(0D0)) :: csensmax ! maximum
    REAL(KIND(0D0)) :: csensavg !mean coverage value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    csens=0. ! csens was already allocated (low dimensional) in bsmatm 

    IF (lip) THEN
       DO j=1,manz
          DO i=1,nanz
             csens(j) = csens(j) + DBLE(sens(i,j)) * &
                  DBLE(sens(i,j)) * wmatd(i)*DBLE(wdfak(i))
          END DO
       END DO
    ELSE IF (ldc) THEN
       DO j=1,manz
          DO i=1,nanz
             csens(j) = csens(j) + sensdc(i,j) * &
                  sensdc(i,j) * wmatd(i)*DBLE(wdfak(i))
          END DO
       END DO
    ELSE
       DO j=1,manz
          DO i=1,nanz
             csens(j) = csens(j) + DCONJG(sens(i,j)) * &
                  sens(i,j) * wmatd(i)*dble(wdfak(i)) 
!!!$ wechselt automatisch wmatdp bei lip
          END DO
       END DO
    ENDIF

!!!$ for normalization
    csensmax = MAXVAL(csens)

    csensavg = SUM (csens) / DBLE(manz)

  END SUBROUTINE bcsens

  subroutine bsmatmreg

!!!$     Unterprogramm belegt die Rauhigkeitsmatrix.
!!!$
!!!$    Andreas Kemna                                        29-Feb-1996
!!!$    Letzte Aenderung   RM                                Jul-2010
!!!$
!!!$....................................................................
!!!$    PROGRAMMINTERNE PARAMETER:
!!!$
!!!$    Variablen zur Beruecksichtigung von Diskontinuitaeten (keine
!!!$    Glaettung in x bzw. z-Richtung)
    INTEGER(KIND = 4) ::   ndis_z,idis_z(3),m,ndis_x,idis_x(4)
    LOGICAL ::     lup,ldown,lleft,lright
    REAL(KIND(0D0)) ::     alfdis

!!!$    Hilfsvariablen
    REAL(KIND(0D0)) ::     dum,dzleft,dzright,xleft,xmean,xright,&
         dxup,dxdown,zup,zmean,zdown

    INTEGER(KIND = 4) ::     i,j,l

!!!$    Hilfsfunction
    INTEGER(KIND = 4) ::     k

    k(i,j) = (i-1) * nx + j

!!!$....................................................................

    IF (.NOT.ALLOCATED (smatm)) ALLOCATE (smatm(manz,3))
    ndis_z = 0
!!!$    ak BAW
!!!$    ak        ndis_z    = 2
    idis_z(1) = 15
    idis_z(2) = 18

    ndis_x = 0
!!!$    ak Bohrloch-Effekt
!!!$    ak        ndis_x    = 4
    idis_x(1) = 4
    idis_x(2) = 6
    idis_x(3) = nx-4
    idis_x(4) = nx-2

!!!$    ak fuer Christoph (Wald)
!!!$    ak        ndis_z    = 3
    idis_z(1) = 6
    idis_z(2) = 12
    idis_z(3) = 19

    IF (.NOT.ALLOCATED (smatm)) ALLOCATE (smatm(manz,3))

!!!$    Rauhigkeitsmatrix auf Null setzen
    do i=1,manz
       do j=1,3
          smatm(i,j) = 0d0
       end do
    end do

    do i=1,nz
       lup   = .true.
       ldown = .true.

       do m=1,ndis_z
          if (i  .eq.idis_z(m)) lup  =.false.
          if (i+1.eq.idis_z(m)) ldown=.false.
       end do

       do j=1,nx
          lleft  = .true.
          lright = .true.

          do m=1,ndis_x
             if (j  .eq.idis_x(m)) lleft =.false.
             if (j+1.eq.idis_x(m)) lright=.false.
          end do

!!!$    Beitrag von Wx^t*Wx zur Rauhigkeitsmatrix
          dzleft  = dabs( sy(snr(nrel(k(i,j),4))) &
               -sy(snr(nrel(k(i,j),1))))
          dzright = dabs( sy(snr(nrel(k(i,j),3))) &
               -sy(snr(nrel(k(i,j),2))))

          xmean = 0d0
          do l=1,4
             xmean = xmean + sx(snr(nrel(k(i,j),l)))
          end do
          xmean = xmean/4d0


          if (j.gt.1) then
             if (lleft) then
                alfdis = 1d0
             else
!!!$    ak
                alfdis = 1d-3
             end if

             xleft = 0d0
             do l=1,4
                xleft = xleft + sx(snr(nrel(k(i,j-1),l)))
             end do

             xleft = xleft/4d0

             smatm(k(i,j),1) = alfdis*alfx * dzleft/dabs(xmean-xleft)
          end if

          if (j.lt.nx) then
             if (lright) then
                alfdis = 1d0
             else
!!!$    ak
                alfdis = 1d-3
             end if

             xright = 0d0
             do l=1,4
                xright = xright + sx(snr(nrel(k(i,j+1),l)))
             end do

             xright = xright/4d0
             dum    = alfdis*alfx * dzright/dabs(xright-xmean)

             smatm(k(i,j),1) = smatm(k(i,j),1) + dum
             smatm(k(i,j),2) = -dum
          end if

!!!$    Beitrag von Wz^t*Wz zur Rauhigkeitsmatrix
          dxup   = dabs( sx(snr(nrel(k(i,j),3))) &
               -sx(snr(nrel(k(i,j),4))))
          dxdown = dabs( sx(snr(nrel(k(i,j),2))) &
               -sx(snr(nrel(k(i,j),1))))

          zmean = 0d0
          do l=1,4
             zmean = zmean + sy(snr(nrel(k(i,j),l)))
          end do
          zmean = zmean/4d0

          if (i.gt.1) then
             if (lup) then
                alfdis = 1d0
             else
                alfdis = 0d0
             end if

             zup = 0d0
             do l=1,4
                zup = zup + sy(snr(nrel(k(i-1,j),l)))
             end do

             zup = zup/4d0
             dum = alfdis*alfz * dxup/dabs(zup-zmean)

             smatm(k(i,j),1) = smatm(k(i,j),1) + dum
          end if

          if (i.lt.nz) then
             if (ldown) then
                alfdis = 1d0
             else
                alfdis = 0d0
             end if

             zdown = 0d0
             do l=1,4
                zdown = zdown + sy(snr(nrel(k(i+1,j),l)))
             end do

             zdown = zdown/4d0
             dum   = alfdis*alfz * dxdown/dabs(zmean-zdown)

             smatm(k(i,j),1) = smatm(k(i,j),1) + dum
             smatm(k(i,j),3) = -dum
          end if

       end do
    end do
  end subroutine bsmatmreg

  SUBROUTINE bsmatmtri      !tri
!!!$
!!!$    Unterprogramm belegt die Rauhigkeitsmatrix....
!!!$    fuer beliebige Triangulierung 
!!!$    Angelehnt an R. Blaschek (2008)
!!!$
!!!$    Copyright by Andreas Kemna 2009
!!!$    
!!!$    Created by Roland Martin                            23-Jun-2009
!!!$    
!!!$    Letzte Aenderung   RM                                    Jul-2010
!!!$
!!!$........................................................................
!!!$    PROGRAMMINTERNE PARAMETER:
!!!$
!!!$    Hilfsvariablen 
    REAL(KIND(0D0)) :: dum    !!!$dummy stores numbers
    INTEGER         :: i,l,k,ik,anz
    REAL(KIND(0D0)) :: edglen !!!$Kantenlaenge
    REAL(KIND(0D0)) :: dist   !!!$Abstand der Schwerpunkte
    REAL(KIND(0D0)) :: sp1(2),sp2(2) !!!$Schwerpunktkoordinaten
    REAL(KIND(0D0)) :: ang    !Winkel fuer anisotrope Glaettung
    REAL(KIND(0D0)) :: alfgeo !Anisotrope (geometrische) Glaettung
!!!$.....................................................................

    IF (.NOT.ALLOCATED (smatm)) ALLOCATE (smatm(manz,smaxs+1),STAT=errnr)
    IF (errnr/=0) THEN
       WRITE (*,'(/a/)')'Allocation problem smatm in bsmatmtri'
       errnr = 97
       RETURN
    END IF

    smatm = 0D0               !!!$initialize smatm

    IF (elanz/=manz)PRINT*,'elanz/=manzSMATMTRI may be wrong'

    DO i=1,elanz

       sp1(1) = espx(i)
       sp1(2) = espy(i)

       DO k=1,smaxs           !!!$jedes flaechenele hat mind einen nachbarn

          ik = MOD(k,smaxs) + 1

          edglen = SQRT((sx(snr(nrel(i,k))) - sx(snr(nrel(i,ik))))**2 + &
               (sy(snr(nrel(i,k))) -  sy(snr(nrel(i,ik))))**2) !!!$edge

          IF (nachbar(i,k)>0) THEN !nachbar existiert 

             sp2(1) = espx(nachbar(i,k)) !!!$schwerpunkt des nachbar elements
             sp2(2) = espy(nachbar(i,k))
!!!$    Geometrischer Teil...

             dist = SQRT((sp1(1) - sp2(1))**2. + (sp1(2) - sp2(2))**2.)

             ang = DATAN2((sp1(2) - sp2(2)),(sp1(1) - sp2(1))) !Winkel

             alfgeo = DSQRT((alfx*DCOS(ang))**2. + (alfz*DSIN(ang))**2.)
             
             dum = edglen / dist * alfgeo
             
             smatm(i,k) = -dum ! set off diagonal of R
             
             smatm(i,smaxs+1) = smatm(i,smaxs+1) + dum ! Main diagonal

          END IF

       END DO
    END DO

  END SUBROUTINE bsmatmtri

  SUBROUTINE bsmatmlma      ! levenberg-marquardt damping
!!!$    
!!!$    Unterprogramm belegt die Glaettungsmtrix
!!!$    (hier Daempfungsmatrix (smatm))
!!!$    
!!!$    Copyright by Andreas Kemna                               2009
!!!$    Erstellt von Roland Martin                               18-Dec-2009
!!!$    Letzte Aenderung RM                                      Jul-2010
!!!$    
!!!$........................................................................
!!!$    PROGRAMMINTERNE PARAMETER:
    REAL(KIND(0D0)) :: csensmax  !Maximale Covarage
    REAL(KIND(0D0)) :: csensavg  !Mittlere Covarage
!!!$     Hilfsvariablen 
    INTEGER         :: i,l,k,j
!!!$.....................................................................

    IF (.NOT.ALLOCATED (smatm)) ALLOCATE (smatm(manz,1),STAT=errnr)
    IF (errnr/=0) THEN
       WRITE (*,'(/a/)')'Allocation problem smatm in bsmatmlma'
       errnr = 97
       RETURN
    END IF

    smatm = 0d0               ! initialize smatm

    IF (ltri==3) THEN

       smatm = 1.0 ! Levenberg Damping

    ELSE

       CALL bcsens (csensmax,csensavg)

       smatm(:,1) = csens

    END IF

  END SUBROUTINE bsmatmlma


  SUBROUTINE bsmatmmgs      !MGS
!!!$    
!!!$  Unterprogramm belegt die Rauhigkeitsmatrix ala 
!!!$ Portniaguine und Zhdanov [1999]
!!!$ Fuer beliebige Triangulierung mit Sensitivitäten gewichtet [Blaschek 2008]
!!!$
!!!$    Copyright by Andreas Kemna 2009
!!!$    
!!!$    Erste Version von Roland Martin                          03-Nov-2009
!!!$    
!!!$    Last edited  RM                                          18-Dec-2009
!!!$    
!!!$........................................................................
!!!$     PROGRAMMINTERNE PARAMETER:
!!!$     Hilfsvariablen 

    REAL(KIND(0D0)) :: dum,dum2
    INTEGER         :: i,j,l,k,ik,anz
    REAL(KIND(0D0)) :: edglen ! Kantenlaenge
    REAL(KIND(0D0)) :: dist ! Abstand der Schwerpunkte
    REAL(KIND(0D0)) :: sp1(2),sp2(2) ! Schwerpunktkoordinaten
    REAL(KIND(0D0)) :: ang    !Winkel fuer anisotrope Glaettung
    REAL(KIND(0D0)) :: csensmax  !Maximale Covarage
    REAL(KIND(0D0)) :: csensavg  !Mittlere Covarage
    REAL(KIND(0D0)) :: alfgeo !Anisotrope Glaettung
    REAL(KIND(0D0)) :: alfmgs !MGS Glaettung
    REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: csens 
!!!$.....................................................................

    errnr = 4
    
    CALL bcsens(csensmax,csensavg)

    IF (csensmax > 1d-12) THEN
       csens = csens / csensmax
    END IF

    WRITE(*,*) 'dum csensavg',dum,csensavg

    IF (.NOT.ALLOCATED(smatm)) ALLOCATE (smatm(manz,smaxs+1),&
         STAT=errnr)

    IF (errnr/=0) THEN
       fetxt = 'Allocation problem WORK in bmcm'
       WRITE (*,'(/a/)')TRIM(fetxt)
       errnr = 97
       RETURN
    END IF

    smatm = 0d0               ! initialize smatm

    DO i=1,elanz ! elanz = flaecheneles

       sp1(1) = espx(i) ! Mittelpunkt des aktuellen Elements
       sp1(2) = espy(i)

       DO k=1,smaxs           ! jedes flaechenele hat mind einen nachbarn

          ik = MOD(k,smaxs) + 1

          edglen = SQRT((sx(snr(nrel(i,k))) - sx(snr(nrel(i,ik))))**2 + &
               (sy(snr(nrel(i,k))) - sy(snr(nrel(i,ik))))**2) 
!!$! edge of i,k and the next..


          IF (nachbar(i,k)>0) THEN !nachbar existiert 

             sp2(1) = espx(nachbar(i,k))
             sp2(2) = espy(nachbar(i,k))
!!!$schwerpunkt des nachbar elements

!!!$    Geometrischer Teil...
             dist = SQRT((sp1(1) - sp2(1))**2. + (sp1(2) - sp2(2))**2.)
             
             ang = DATAN2((sp1(2) - sp2(2)),(sp1(1) - sp2(1))) !neu
             
             alfgeo = DSQRT((alfx*DCOS(ang))**2. + (alfz*DSIN(ang))**2.)
!!!$    MGS Teil
             dum = CDABS(sigma(i) - sigma(nachbar(i,k))) / dist
!!!$    Modell nur im zaehler nicht im Nenner -> fred fragen

!!!$    \int \frac{(\nabla m_{ij})^2}{(\nabla m_{ij})^2+\beta^2}\;dA
!!!$    -> (m_i-m_{i+1})^2 \frac{\Delta z_i}{\Delta x_i}
!!!$    wobei (m_i-m_{i+1})^2 rausgezogen wird und spaeter 
!!!$    als Matrix Vektor Produkt berechnet wird
!!!$    (m_i-m_{i+1})^2 \frac{\Delta z_i}{\Delta x_i} 
!!!$    -> smatm(i) = \frac{\Delta z_i}{\Delta x_i} * geometrischem Teil 
!!!$    von anisotroper Wichtung
             IF (ltri == 5) THEN !!!$reines MGS

                dum = dum**2. + betamgs**2.
                dum = alfgeo * edglen / dist / dum

             ELSE IF (ltri == 6) THEN !!!$sensitivitaetswichtung 1 von RM
!!!$    f(i,k) = 1 + g(i) + g(k)
                dum2 = 1d0 + DABS(DLOG10(csens(i))) + &
                     DABS(DLOG10(csens(nachbar(i,k))))
!!!$    dum2 = f(i,k)^2
                dum2 = dum2**2.
!!!$    dum = grad(m)^2 + (\beta/f(i,k)^2)^2
                dum = dum**2. + (betamgs / dum2)**2.
!!!$    dum = \alpha_{xz} * \Delta z / \Delta x / f(i,k)^2 / 
!!!$    grad(m)^2 + (\beta/f(i,k)^2)^2
                dum = alfgeo * edglen / dist /dum2 / dum

             ELSE IF (ltri == 7) THEN !!!$sensitivitaetswichtung 2 von RM

!!!$    f(i,k) = 1 + (g(i) + g(k))/mean(g)
                dum2 = 1d0 + DABS((DLOG10(csens(i))) + &
                     DABS(DLOG10(csens(nachbar(i,k))))) / csensavg
!!!$    dum2 = f(i,k)^2
                dum2 = dum2**2.
!!!$    dum = grad(m)^2 + (\beta/f(i,k)^2)^2
                dum = dum**2. + (betamgs / dum2)**2.
!!!$    dum = \alpha_{xz} * \Delta z / \Delta x / f(i,k)^2 / 
!!!$    grad(m)^2 + (\beta/f(i,k)^2)^2
                dum = alfgeo * edglen / dist /dum2 / dum

             ELSE IF (ltri == 8) THEN !!!$sensitivitaetswichtung von RB
                
!!!$    der folgende code wurde mir so ueberliefert... 
!!!$ kam von RB aber keine ahnung was das genau macht
                dum = dum * (1d0 + 0.2d0 * (DABS( DLOG10(csens(i)) + & 
                     DLOG10(csens(nachbar(i,k))) ) ))
                
                alfmgs = 1d0 - dum**2. / (dum**2. + betamgs**2.)
                dum =  edglen * alfgeo * alfmgs

             ELSE IF (ltri == 9) THEN

                dum = dum * (1d0 + 0.2d0 * (DABS( DLOG10(csens(i)) + &
                     DLOG10(csens(nachbar(i,k))) ) / csensavg ))

                alfmgs = 1d0 - dum**2. / (dum**2. + betamgs**2.)
                dum =  edglen * alfgeo * alfmgs

             END IF

!!!$    nun glaettung belegen
             smatm(i,k) = -dum !!!$neben Diagonale
             smatm(i,smaxs+1) = smatm(i,smaxs+1) + dum !Hauptdiagonale

          END IF
       END DO
    END DO

    DEALLOCATE (csens)

    errnr = 0
999 RETURN

  END SUBROUTINE bsmatmmgs

  SUBROUTINE bsmatmtv      !betatv
!!!$    
!!!$    Unterprogramm belegt die Rauhigkeitsmatrix mit total variance
!!!$    fuer beliebige Triangulierung 
!!!$
!!!$    Copyright by Andreas Kemna 2009
!!!$    
!!!$    Created by Roland Martin                          23-Nov-2009
!!!$    
!!!$    Letzte Aenderung   RM                                    23-Nov-2009
!!!$    
!!!$........................................................................
!!!$   PROGRAMMINTERNE PARAMETER:
!!!$   Hilfsvariablen 
    REAL(KIND(0D0)) :: dum
    INTEGER         :: i,j,l,k,ik,anz
    REAL(KIND(0D0)) :: edglen !!!$Kantenlaenge
    REAL(KIND(0D0)) :: dist   !!!$Abstand der Schwerpunkte
    REAL(KIND(0D0)) :: sp1(2),sp2(2) !!!$Schwerpunktkoordinaten
    REAL(KIND(0D0)) :: ang    !Winkel fuer anisotrope Glaettung
    REAL(KIND(0D0)) :: alfgeo !Anisotrope (geometrische) Glaettung
    REAL(KIND(0D0)) :: alftv  !TV Glaettung
!!!$.....................................................................

    
    IF (.NOT.ALLOCATED(smatm)) ALLOCATE (smatm(manz,smaxs+1))
    smatm = 0d0               !!!$initialize smatm

    DO i=1,elanz
       sp1(1) = espx(i) !!!$Mittelpunkt des aktuellen Elements
       sp1(2) = espy(i)

       DO k=1,smaxs           !!!$jedes flaechenele hat mind einen nachbarn

          ik = MOD(k,smaxs) + 1

          edglen = SQRT((sx(snr(nrel(i,k))) - sx(snr(nrel(i,ik))))**2 + &
               (sy(snr(nrel(i,k))) - sy(snr(nrel(i,ik))))**2) !!!$edge


          IF (nachbar(i,k)>0) THEN !nachbar existiert 

             sp2(1) = espx(nachbar(i,k)) !!!$schwerpunkt des nachbar elements
             sp2(2) = espy(nachbar(i,k))

!!!$   Geometrischer Teil...
             dist = SQRT((sp1(1) - sp2(1))**2. + (sp1(2) - sp2(2))**2.)

             ang = DATAN2((sp1(2) - sp2(2)),(sp1(1) - sp2(1))) !neu

             alfgeo = DSQRT((alfx*DCOS(ang))**2. + (alfz*DSIN(ang))**2.)

             alftv = edglen / dist * alfgeo

!!!$   Total variance
             dum = SQRT(alftv**2. + betamgs**2.)
!!!$   nun glaettung belegen

             smatm(i,k) = -dum !!!$ off diagonal

             smatm(i,smaxs+1) = smatm(i,smaxs+1) + dum !main diagonal

          END IF

       END DO
    END DO

  END SUBROUTINE bsmatmtv

  SUBROUTINE bsmatmsto
!!!$
!!!$    Unterprogramm belegt die Kovarianzmatrix.   
!!!$    Neue Regularisierungsmatrix (stoch. Kovarianzmatrix).
!!!$
!!!$    Copyright by Andreas Kemna                              2009
!!!$    Created by Anastasia August / Roland Martin             03-Apr-2009
!!!$
!!!$    Last edited RM                                          Jul-2010
!!!$....................................................................
!!!$    Hilfsmatrix
    REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: work
!!!$    Schwerpunktskoordinaten der Flaechenelemente ij
    REAL(KIND(0D0)) :: h,sd_el
!!!$    Korrelation lengths and variance (var)
    REAL(KIND(0D0))      :: hx,hy,var
!!!$    gibt es evtl schon eine inverse?
    logical              :: ex,exc        
!!!$    Hilfsvariablen
    integer              :: i,j,l,ifp,c1
!!!$    smatm file name
    CHARACTER(124)        :: fsmat
!!!$    clearscreen
    CHARACTER(80)        :: csz


    errnr = 1
    CALL get_unit(ifp)

    var = 1.

    fsmat = ramd(1:lnramd)//slash(1:1)//'inv.smatmi'

    DO i=1,79
       csz(i:i+1)=' '
    END DO

    WRITE (*,'(A80)')ACHAR(13)//TRIM(csz)

    WRITE (*,'(A,1X,F6.2,1X,A)')ACHAR(13)//'Speicher fuer model '//&
         'covariance: ',REAL ((manz**2*8.)/(1024.**3)),'GB'
    
    IF (.NOT.ALLOCATED (smatm))ALLOCATE (smatm(manz,manz),STAT=errnr)
    IF (errnr/=0) THEN
       WRITE (*,'(/a/)')'Allocation problem smatm in bsmatmsto'
       errnr = 97
       RETURN
    END IF

!!!$    Belege die Matrix

    smatm = var

    INQUIRE(FILE=fsmat,EXIST=ex) !!!$already an inverse c_m ?

    IF (ex) THEN

       WRITE (*,'(a)',ADVANCE='no')'checking '//fsmat
       OPEN (ifp,FILE=fsmat,STATUS='old',ACCESS='sequential',&
            FORM='unformatted')
       READ (ifp) i
       IF (i == manz) THEN
          WRITE(*,'(a)')'ok!'
          READ (ifp) smatm
       END IF
       CLOSE (ifp)

       errnr = 0

    ELSE

       DO i = 1 , manz
          WRITE (*,'(a,1X,F6.2,A)',ADVANCE='no')ACHAR(13)//'cov/',&
               REAL(i*(100./manz)),'%'

          DO j = i+1 , manz   !!!$fills upper triangle

             hx = (espx(i) - espx(j))
             hy = (espy(i) - espy(j))

             smatm(i,j) = mcova(hx,hy,var)

             smatm(j,i) = smatm(i,j) !!!$upper triangle

          END DO
       END DO


       PRINT*,'bestimme nun C_m^-1'
!!!$    Berechne nun die Inverse der Covarianzmatrix!!!
       IF (lgauss) THEN
          PRINT*,'   Gauss elemination ... '
          CALL gauss_dble(smatm,manz,errnr)
          IF (errnr/=0) THEN
             fetxt='there was something wrong..'
             PRINT*,'Zeile(',abs(errnr),')::',smatm(abs(errnr),:)
             errnr = 108
             RETURN
          END IF
       ELSE                   !!!$default..
          WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//'Factorization...'
          ALLOCATE (work(manz))
          CALL CHOLD(smatm,work,manz,errnr)
          IF (errnr/=0) THEN
             fetxt='CHOLD smatm :: matrix not pos definite..'
             PRINT*,'Zeile(',abs(errnr),')'
             errnr = 108
             RETURN
          END IF
          WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//'Inverting...'
          CALL LINVD(smatm,work,manz)
          DEALLOCATE (work)
          WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//'Filling upper C_m...'
          DO i= 1, manz
             WRITE (*,'(A,1X,F6.2,A)',ADVANCE='no')ACHAR(13)//ACHAR(9)//&
                  ACHAR(9)//ACHAR(9)//'/ ',REAL( i * (100./manz)),'%'
             DO j = 1, i
                smatm(i,j) = smatm(j,i)
             END DO
          END DO
       END IF

       IF (errnr==0) THEN
          WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//'got inverse and write out'
          OPEN (ifp,FILE=fsmat,STATUS='replace',ACCESS='sequential',&
               FORM='unformatted')
          WRITE (ifp) manz
          WRITE (ifp) smatm
          CLOSE (ifp)
       ELSE
          PRINT*,'got NO inverse'
          errnr = 108
          RETURN
       END IF
    END IF

  END SUBROUTINE bsmatmsto

END MODULE bsmatm_mod

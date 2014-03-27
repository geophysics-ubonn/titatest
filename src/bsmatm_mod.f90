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
  USE alloci , ONLY : sens,sensdc,smatm,prec
  USE femmod , ONLY : fak,ldc
  USE elemmod, ONLY : smaxs,sx,sy,espx,espy,nrel,snr,elanz,nachbar
  USE invmod , ONLY : lfpi,par,wmatd,wdfak
  USE errmod , ONLY : errnr,fetxt
  USE konvmod , ONLY : ltri,lgauss,lam,nx,nz,alfx,alfz,betamgs,lverb,lverb_dat
  USE modelmod , ONLY : manz
  USE datmod , ONLY : nanz
  USE errmod, ONLY : errnr,fetxt
  USE sigmamod , ONLY : sigma
  USE ompmod
  USE variomodel 
  USE pathmod
  IMPLICIT NONE

  REAL(prec),DIMENSION(:),ALLOCATABLE,PRIVATE :: csens 

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

  SUBROUTINE bsmatm(it,l_bsmat)
!!!$
!!!$ This sub is the control unit of the smatm calculation
!!!$
    INTEGER (KIND = 4 ),INTENT(IN) :: it
    INTEGER (KIND = 4 )   :: c1
    LOGICAL            ,INTENT (INOUT) :: l_bsmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    errnr = 2

!!!$ pre query if smatm is  ever to be recalculated again..

    IF (ltri > 4 .AND. ltri < 10) THEN
       
       IF (betamgs < 0d0) THEN
          PRINT*
          PRINT*,'Fixing Regularization matrix'
          PRINT*
          betamgs = ABS(betamgs)
          l_bsmat = .FALSE.
       ELSE
          l_bsmat = .TRUE.
       END IF
    ELSE IF (ltri == 4) THEN ! for levenberg-marquardt JTJ scaling
       l_bsmat = .TRUE.
    ELSE
       l_bsmat = .FALSE. ! calculate just once
    END IF

!!!$    get time
    CALL TIC(c1)

    IF (.NOT.ALLOCATED(csens)) ALLOCATE (csens(manz),STAT=errnr)
    IF (errnr/=0) THEN
       fetxt = 'Allocation problem csens in bsmatmmgs'
       WRITE (*,'(/a/)')TRIM(fetxt)
       errnr = 97
       RETURN
    END IF
    
    IF (ltri == 0) THEN
       
       WRITE (*,'(a)')' rectangular smooth regularization'
       CALL bsmatmreg
       
    ELSE IF (ltri == 1.OR.ltri == 2) THEN
       
       WRITE (*,'(a)')' triangular smooth regularization'
       CALL bsmatmtri
       
    ELSE IF (ltri == 3.OR.ltri==4) THEN
       
       IF (it == 1) THEN
          IF (ltri == 3) WRITE (*,'(a)') &
               ' Levenberg damping (LA)'
          IF (ltri == 4) WRITE (*,'(a)') &
               ' Marquardt-Levenberg damping (LMA)'
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

    fetxt = 'C_m calculation time'
    CALL TOC(c1,fetxt)

  END SUBROUTINE bsmatm

!!! 
  SUBROUTINE bcsens (csensmax,csensavg)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!$ This subroutine calculates the squared coverage or diag{A^TC_d^{-1}A}!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER :: i,j
    REAL(prec) :: csensmax ! maximum
    REAL(prec) :: csensavg !mean coverage value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    csens=0d0 ! csens was already allocated (low dimensional) in bsmatm 

    !$OMP PARALLEL DEFAULT (none) &
    !$OMP SHARED (csens,wmatd,sens,sensdc,wdfak,lfpi,ldc,manz,nanz) &
    !$OMP PRIVATE(j,i)
    !$OMP DO SCHEDULE (GUIDED,CHUNK_0)
    DO j=1,manz
       DO i=1,nanz
          IF (lfpi) THEN
             csens(j) = csens(j) + REAL(sens(i,j)) * &
                  REAL(sens(i,j)) * wmatd(i)*REAL(wdfak(i))
          ELSE IF (ldc) THEN
             csens(j) = csens(j) + sensdc(i,j) * &
                  sensdc(i,j) * wmatd(i)*REAL(wdfak(i))
!!!$ wechselt automatisch wmatdp bei lfpi
          ELSE
             csens(j) = csens(j) + CONJG(sens(i,j)) * &
                  sens(i,j) * wmatd(i)*REAL(wdfak(i)) 
          ENDIF
       END DO
    END DO
    !$OMP END PARALLEL

!!!$ for normalization
    csensmax = MAXVAL(csens)

    csensavg = SUM (csens) / REAL(manz)

  END SUBROUTINE bcsens

  SUBROUTINE bsmatmreg

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
    REAL(prec) ::     alfdis

!!!$    Hilfsvariablen
    REAL(prec) ::     dum,dzleft,dzright,xleft,xmean,xright,&
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
    smatm = 0d0

    DO i=1,nz
       lup   = .TRUE.
       ldown = .TRUE.

       DO m=1,ndis_z
          IF (i  .EQ.idis_z(m)) lup  =.FALSE.
          IF (i+1.EQ.idis_z(m)) ldown=.FALSE.
       END DO

       DO j=1,nx
          lleft  = .TRUE.
          lright = .TRUE.

          DO m=1,ndis_x
             IF (j  .EQ.idis_x(m)) lleft =.FALSE.
             IF (j+1.EQ.idis_x(m)) lright=.FALSE.
          END DO

!!!$    Beitrag von Wx^t*Wx zur Rauhigkeitsmatrix
          dzleft  = ABS( sy(snr(nrel(k(i,j),4))) &
               -sy(snr(nrel(k(i,j),1))))
          dzright = ABS( sy(snr(nrel(k(i,j),3))) &
               -sy(snr(nrel(k(i,j),2))))

          xmean = 0d0
          DO l=1,4
             xmean = xmean + sx(snr(nrel(k(i,j),l)))
          END DO
          xmean = xmean/4d0


          IF (j.GT.1) THEN
             IF (lleft) THEN
                alfdis = 1d0
             ELSE
!!!$    ak
                alfdis = 1d-3
             END IF

             xleft = 0d0
             DO l=1,4
                xleft = xleft + sx(snr(nrel(k(i,j-1),l)))
             END DO

             xleft = xleft/4d0

             smatm(k(i,j),1) = alfdis*alfx * dzleft/ABS(xmean-xleft)
          END IF

          IF (j.LT.nx) THEN
             IF (lright) THEN
                alfdis = 1d0
             ELSE
!!!$    ak
                alfdis = 1d-3
             END IF

             xright = 0d0
             DO l=1,4
                xright = xright + sx(snr(nrel(k(i,j+1),l)))
             END DO

             xright = xright/4d0
             dum    = alfdis*alfx * dzright/ABS(xright-xmean)

             smatm(k(i,j),1) = smatm(k(i,j),1) + dum
             smatm(k(i,j),2) = -dum
          END IF

!!!$    Beitrag von Wz^t*Wz zur Rauhigkeitsmatrix
          dxup   = ABS( sx(snr(nrel(k(i,j),3))) &
               -sx(snr(nrel(k(i,j),4))))
          dxdown = ABS( sx(snr(nrel(k(i,j),2))) &
               -sx(snr(nrel(k(i,j),1))))

          zmean = 0d0
          DO l=1,4
             zmean = zmean + sy(snr(nrel(k(i,j),l)))
          END DO
          zmean = zmean/4d0

          IF (i.GT.1) THEN
             IF (lup) THEN
                alfdis = 1d0
             ELSE
                alfdis = 0d0
             END IF

             zup = 0d0
             DO l=1,4
                zup = zup + sy(snr(nrel(k(i-1,j),l)))
             END DO

             zup = zup/4d0
             dum = alfdis*alfz * dxup/ABS(zup-zmean)

             smatm(k(i,j),1) = smatm(k(i,j),1) + dum
          END IF

          IF (i.LT.nz) THEN
             IF (ldown) THEN
                alfdis = 1d0
             ELSE
                alfdis = 0d0
             END IF

             zdown = 0d0
             DO l=1,4
                zdown = zdown + sy(snr(nrel(k(i+1,j),l)))
             END DO

             zdown = zdown/4d0
             dum   = alfdis*alfz * dxdown/ABS(zmean-zdown)

             smatm(k(i,j),1) = smatm(k(i,j),1) + dum
             smatm(k(i,j),3) = -dum
          END IF

       END DO
    END DO
  END SUBROUTINE bsmatmreg

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
    REAL(prec) :: dum    !!!$dummy stores numbers
    INTEGER         :: i,k,ik
    REAL(prec) :: edglen !!!$Kantenlaenge
    REAL(prec) :: dist   !!!$Abstand der Schwerpunkte
    REAL(prec) :: sp1(2),sp2(2) !!!$Schwerpunktkoordinaten
    REAL(prec) :: ang    !Winkel fuer anisotrope Glaettung
    REAL(prec) :: alfgeo !Anisotrope (geometrische) Glaettung
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

          ik = MOD(k,smaxs) + 1 !!! associates the next node, or itself

          edglen = SQRT((sx(snr(nrel(i,k))) - sx(snr(nrel(i,ik))))**2d0 + &
               (sy(snr(nrel(i,k))) -  sy(snr(nrel(i,ik))))**2d0) !!!$edge

          IF (nachbar(i,k)>0) THEN !nachbar an der Kante existiert 

             sp2(1) = espx(nachbar(i,k)) !!!$ center point of element
             sp2(2) = espy(nachbar(i,k))
!!!$    Geometrical part...

             dist = SQRT((sp1(1) - sp2(1))**2d0 + (sp1(2) - sp2(2))**2d0)

             ang = ATAN2((sp1(2) - sp2(2)),(sp1(1) - sp2(1))) !Angle

             alfgeo = SQRT((alfx*COS(ang))**2d0 + (alfz*SIN(ang))**2d0)
             
             dum = edglen / dist * alfgeo ! proportional contribution of integrated cell
             
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
    REAL(prec) :: csensmax  !Maximale Covarage
    REAL(prec) :: csensavg  !Mittlere Covarage
    INTEGER :: j
!!!$.....................................................................

    IF (.NOT.ALLOCATED (smatm)) ALLOCATE (smatm(manz,1),STAT=errnr)
    IF (errnr/=0) THEN
       WRITE (*,'(/a/)')'Allocation problem smatm in bsmatmlma'
       errnr = 97
       RETURN
    END IF

    smatm = 0d0               ! initialize smatm

    IF (ltri==3) THEN

       smatm = 1d0 ! Levenberg Damping

    ELSE

       CALL bcsens (csensmax,csensavg)

       DO j = 1,manz

          smatm(j,1) = csensmax/csens(j)

       END DO

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

    REAL(prec) :: dum,dum2 ! helpers
    REAL(prec) :: mgrad,sqmgrad ! model gradient and squared model grad
    INTEGER         :: i,k,ik
    REAL(prec) :: edglen ! Kantenlaenge
    REAL(prec) :: dist ! Abstand der Schwerpunkte
    REAL(prec) :: sp1(2),sp2(2) ! Schwerpunktkoordinaten
    REAL(prec) :: ang    !Winkel fuer anisotrope Glaettung
    REAL(prec) :: csensmax  !Maximale Covarage
    REAL(prec) :: csensavg  !Mittlere Covarage
    REAL(prec) :: alfgeo !Anisotrope Glaettung
    REAL(prec) :: alfmgs !MGS Glaettung
!!!$.....................................................................

    errnr = 4
    
    CALL bcsens(csensmax,csensavg)

    IF (csensmax > 1d-12) THEN
       csens = csens / csensmax
    END IF

    PRINT*,'csensavg/csensmax',csensavg,'/',csensmax

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

          ik = MOD(k,smaxs) + 1 !!! associates the next node, or itself

          edglen = SQRT((sx(snr(nrel(i,k))) - sx(snr(nrel(i,ik))))**2d0 + &
               (sy(snr(nrel(i,k))) - sy(snr(nrel(i,ik))))**2d0) 
!!$! edge of i,k and the next..

          IF (nachbar(i,k)>0) THEN !nachbar existiert 

!!!$schwerpunkt des nachbar elements
             sp2(1) = espx(nachbar(i,k))!!!$ center point of element
             sp2(2) = espy(nachbar(i,k))

!!!$    Geometrical part...
             ! distance of the mid points
             dist = SQRT((sp1(1) - sp2(1))**2d0 + (sp1(2) - sp2(2))**2d0)
!!$! including anisotropy!
!angle to horizon
             ang = ATAN2((sp1(2) - sp2(2)),(sp1(1) - sp2(1)))
!!!$ geometrical contribution... (as smooth regularization..)
! projected effective contribution due to anisotropic regu
             alfgeo = SQRT((alfx*COS(ang))**2d0 + (alfz*SIN(ang))**2d0)

!!!$ Model value gradient (\nabla m)

!!! TODO
             mgrad = ABS(sigma(i) - sigma(nachbar(i,k))) / dist
             sqmgrad = mgrad * mgrad
!!!$ TODO
!!!$    MGS Teil
!!!$
!!!$    \int \frac{(\nabla m_{ij})^2}{(\nabla m_{ij})^2+\beta^2}\;dA
!!!$    -> (m_i-m_{i+1})^2 \frac{\Delta z_i}{\Delta x_i}
!!!$            !!!! ATTENTION !!!!
!!!$ The squared model gradient in the 
!!!$ nominator of the stabilizer,  i.e. (m_i-m_{i+1})^2 
!!!$ !!!   IS EVALUATED LATER ON AS MATRIX VECTOR PRODUCT   !!!
!!!$ for now we have to deal with the denominator stuff only at this point!!
!!!$ The gemoetrical part is than reduced to
!!!$ \frac{\Delta z_i}{\Delta x_i} which is edglen / dist!!!
!!!$ -> smatm(i) = \frac{\Delta z_i}{\Delta x_i} * geometrical part 
!!!$  of anisotropy
             IF (ltri == 5) THEN !!!$reines MGS

                dum = sqmgrad + betamgs**2d0
! proportional contribution of integrated cell
                dum = alfgeo * edglen / dist / dum

             ELSE IF (ltri == 6) THEN !!!$sensitivitaetswichtung 1 von RM
!!!$    f(i,k) = 1 + g(i) + g(k)
                dum2 = 1d0 + ABS(LOG10(csens(i))) + &
                     ABS(LOG10(csens(nachbar(i,k))))
!!!$    dum2 = f(i,k)^2
                dum2 = dum2**2d0
!!!$    dum = grad(m)^2 + (\beta/f(i,k)^2)^2
                dum = sqmgrad + (betamgs / dum2)**2d0
!!!$    dum = \alpha_{xz} * \Delta z / \Delta x / f(i,k)^2 / 
!!!$    grad(m)^2 + (\beta/f(i,k)^2)^2
                dum = alfgeo * edglen / dist / dum2 / dum

             ELSE IF (ltri == 7) THEN !!!$sensitivitaetswichtung 2 von RM

!!!$    f(i,k) = 1 + (g(i) + g(k))/mean(g)
                dum2 = 1d0 + ABS((LOG10(csens(i))) + &
                     ABS(LOG10(csens(nachbar(i,k))))) / csensavg
!!!$    dum2 = f(i,k)^2
                dum2 = dum2**2d0
!!!$    dum = grad(m)^2 + (\beta/f(i,k)^2)^2
                dum = sqmgrad + (betamgs / dum2)**2d0
!!!$    dum = \alpha_{xz} * \Delta z / \Delta x / f(i,k)^2 / 
!!!$    grad(m)^2 + (\beta/f(i,k)^2)^2
                dum = alfgeo * edglen / dist / dum2 / dum

             ELSE IF (ltri == 8) THEN !!!$sensitivitaetswichtung von RB
                
!!!$    der folgende code wurde mir so ueberliefert... 
!!!$ kam von RB aber keine ahnung was das genau macht
                dum = mgrad * (1d0 + 0.2d0 * (ABS( LOG10(csens(i)) + & 
                     LOG10(csens(nachbar(i,k))) ) ))
                
                alfmgs = 1d0 - dum**2d0 / (dum**2d0 + betamgs**2d0)
                dum =  edglen * alfgeo * alfmgs

             ELSE IF (ltri == 9) THEN

                dum = mgrad * (1d0 + 0.2d0 * (ABS( LOG10(csens(i)) + &
                     LOG10(csens(nachbar(i,k))) ) / csensavg ))

                alfmgs = 1d0 - dum**2d0 / (dum**2d0 + betamgs**2d0)
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
    REAL(prec) :: dum
    INTEGER         :: i,k,ik
    REAL(prec) :: edglen !!!$Kantenlaenge
    REAL(prec) :: dist   !!!$Abstand der Schwerpunkte
    REAL(prec) :: sp1(2),sp2(2) !!!$Schwerpunktkoordinaten
    REAL(prec) :: ang    !Winkel fuer anisotrope Glaettung
    REAL(prec) :: alfgeo !Anisotrope (geometrische) Glaettung
    REAL(prec) :: alftv  !TV Glaettung
!!!$.....................................................................

    
    IF (.NOT.ALLOCATED(smatm)) ALLOCATE (smatm(manz,smaxs+1))
    smatm = 0d0               !!!$initialize smatm

    DO i=1,elanz
       sp1(1) = espx(i) !!!$Mittelpunkt des aktuellen Elements
       sp1(2) = espy(i)

       DO k=1,smaxs           !!!$jedes flaechenele hat mind einen nachbarn

          ik = MOD(k,smaxs) + 1

          edglen = SQRT((sx(snr(nrel(i,k))) - sx(snr(nrel(i,ik))))**2d0 + &
               (sy(snr(nrel(i,k))) - sy(snr(nrel(i,ik))))**2d0) !!!$edge


          IF (nachbar(i,k)>0) THEN !nachbar existiert 

             sp2(1) = espx(nachbar(i,k)) !!!$schwerpunkt des nachbar elements
             sp2(2) = espy(nachbar(i,k))

!!!$   Geometrischer Teil...
             dist = SQRT((sp1(1) - sp2(1))**2d0 + (sp1(2) - sp2(2))**2d0)

             ang = ATAN2((sp1(2) - sp2(2)),(sp1(1) - sp2(1))) !neu

             alfgeo = SQRT((alfx*COS(ang))**2d0 + (alfz*SIN(ang))**2d0)

             alftv = edglen / dist * alfgeo

!!!$   Total variance
             dum = SQRT(alftv**2d0 + betamgs**2d0)
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
    REAL(prec),DIMENSION(:),ALLOCATABLE :: work
    REAL(prec),DIMENSION(:,:),ALLOCATABLE :: myold,proof
!!!$    Korrelation lengths, variance (var) and nugget
    REAL(prec)      :: hx,hy,var,nugget
    REAL                 :: epsi
!!!$    gibt es evtl schon eine inverse?
    LOGICAL              :: ex
!!!$    Hilfsvariablen
    INTEGER              :: i,j,ifp
!!!$    smatm file name
    CHARACTER(124)        :: fsmat
!!!$    clearscreen
    CHARACTER(80)        :: csz

    epsi=EPSILON(epsi)
    IF (lverb) WRITE (*,*)'Epsilon smatm::',epsi

    errnr = 1
    CALL get_unit(ifp)

    var = 1d0
    nugget = 1d-4

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
    IF (.NOT.ALLOCATED (proof).AND.lverb) THEN
       ALLOCATE (proof(manz,manz),STAT=errnr)
       IF (errnr/=0) THEN
          WRITE (*,'(/a/)')'Allocation problem smatm in bsmatmsto'
          errnr = 97
          RETURN
       END IF
       ALLOCATE (myold(manz,manz),STAT=errnr)
       IF (errnr/=0) THEN
          WRITE (*,'(/a/)')'Allocation problem smatm in bsmatmsto'
          errnr = 97
          RETURN
       END IF
    END IF
!!!$    Belege die Matrix

    smatm = 0d0

    INQUIRE(FILE=fsmat,EXIST=ex) !!!$already an inverse c_m ?

    IF (ex) THEN

       WRITE (*,'(a)',ADVANCE='no')'checking '//fsmat
       OPEN (ifp,FILE=fsmat,STATUS='old',ACCESS='sequential',&
            FORM='unformatted')
       READ (ifp) i
       IF (i == manz) THEN
          WRITE(*,'(t40,a)')'ok!'
          READ (ifp) smatm
       END IF
       CLOSE (ifp)

       errnr = 0

    ELSE

       !$OMP PARALLEL DEFAULT (none) &
       !$OMP SHARED (smatm,manz,epsi,lverb,ifp,espx,espy,var) &
       !$OMP PRIVATE (i,j,hx,hy)
       !$OMP DO SCHEDULE (GUIDED,CHUNK_0)
       DO i = 1 , manz
          IF (lverb) WRITE (*,'(a,t25,F6.2,A,t70,a)',ADVANCE='no')ACHAR(13)//&
               'cov/',REAL(i*(100./manz)),'%',''

          smatm(i,i) = var ! nugget (=variance) effect on the main
!!!$ R(h) = C_0 \delta(h)  = 
!!!$ \begin{case} C_0 & \mbox{if}\;h = 0 \\ 0 & else \end{case}

          DO j = i+1 , manz   !!!$fills upper triangle

             hx = (espx(i) - espx(j)) !main point differences
             hy = (espy(i) - espy(j))

             smatm(i,j) = mcova(hx,hy,var) ! compute covariance

             smatm(j,i) = smatm(i,j) ! lower triangle

          END DO
       END DO
       !$OMP END PARALLEL

       IF (lverb_dat) THEN
          fetxt = 'cm0.dat'
          PRINT*,'writing '//TRIM(fetxt)
          OPEN (ifp,FILE=TRIM(fetxt),STATUS='replace',&
               ACCESS='sequential',FORM='formatted')
          DO i = 1,manz
             WRITE (ifp,*)espx(i),espy(i),(smatm(i,j),j=i,manz)
          END DO
          CLOSE (ifp)
       END IF

       IF (lverb) myold = smatm

!!!$    Berechne nun die Inverse der Covarianzmatrix!!!
       IF (lgauss) THEN
          PRINT*,'   Gauss elemination ... '
          CALL gauss_REAL(smatm,manz,errnr)
          IF (errnr/=0) THEN
             fetxt='there was something wrong..'
             PRINT*,'Zeile(',ABS(errnr),')::',smatm(ABS(errnr),:)
             errnr = 108
             RETURN
          END IF
       ELSE                   !!!$default..
          WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//'Factorization...'
          ALLOCATE (work(manz))
          CALL CHOLD(smatm,work,manz,errnr,lverb)
          IF (errnr/=0) THEN
             fetxt='CHOLD smatm :: matrix not pos definite..'
             PRINT*,'Zeile(',ABS(errnr),')'
             errnr = 108
             RETURN
          END IF
          WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//'Inverting...'
          CALL LINVD(smatm,work,manz,lverb)
          DEALLOCATE (work)
          !$OMP PARALLEL DEFAULT (none) &
          !$OMP SHARED (smatm,manz,epsi,lverb) &
          !$OMP PRIVATE (i,j)
          !$OMP DO SCHEDULE (GUIDED,CHUNK_0)
          DO i= 1, manz
             IF (lverb) WRITE (*,'(A,t25,F6.2,A)',ADVANCE='no')ACHAR(13)//&
                  'Filling upper C_m',REAL( i * (100./manz)),'%'
             DO j = 1, i - 1

                smatm(i,j) = smatm(j,i)

             END DO
          END DO
          !$OMP END PARALLEL
       END IF

       IF (lverb) THEN

          !$OMP WORKSHARE
          proof = MATMUL(smatm,myold)
          !$OMP END WORKSHARE

          DEALLOCATE (myold)
          DO i=1,manz
             IF (ABS(proof(i,i) - 1d0) > 0.1) PRINT*,'bad approximation at parameter'&
                  ,i,proof(i,i)
          END DO

       END IF

       IF (errnr == 0) THEN
          WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//'got inverse'
          IF (lverb) THEN
             WRITE (*,'(a)',ADVANCE='no')' .. write out'//TRIM(fsmat)
             OPEN (ifp,FILE=fsmat,STATUS='replace',ACCESS='sequential',&
                  FORM='unformatted')
             WRITE (ifp) manz
             WRITE (ifp) smatm
             CLOSE (ifp)
          END IF

       ELSE

          PRINT*,'got NO inverse'
          errnr = 108
          RETURN

       END IF

!!!!$ verbose output of inverse smatm data
       IF (lverb_dat) THEN
          fetxt = 'cm0_inv.dat'
          PRINT*,'writing '//TRIM(fetxt)
          OPEN (ifp,FILE=TRIM(fetxt),STATUS='replace',&
               ACCESS='sequential',FORM='formatted')
          DO i = 1,manz
             WRITE (ifp,*)espx(i),espy(i),(smatm(i,j),j=i,manz)
          END DO
          CLOSE (ifp)
       END IF


    END IF
    IF (ALLOCATED(proof)) DEALLOCATE (proof)

  END SUBROUTINE bsmatmsto

END MODULE bsmatm_mod

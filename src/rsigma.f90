SUBROUTINE rsigma(kanal,datei)

!!!$     Unterprogramm zum Einlesen der Widerstandsverteilung aus 'datei'.

!!!$     Andreas Kemna                                            20-Dec-1993
!!!$     Letzte Aenderung   07-Nov-1997

!!!$.....................................................................

  USE make_noise
  USE alloci,ONLY:rnd_r,rnd_p
  USE femmod
  USE invmod
  USE sigmamod
  USE modelmod
  USE elemmod
  USE errmod
  USE konvmod

  IMPLICIT NONE

!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Kanalnummer
  INTEGER (KIND=4) ::    kanal

!!!$     Datei
  CHARACTER (80)   ::    datei

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsvariablen
  INTEGER (KIND=4) ::     i,idum,ifp1,ifp2
  REAL(KIND(0D0))  ::     bet,pha,eps_r,eps_p,bet_ref,pha_ref
!!!$ Pi
  REAL (KIND(0D0)) ::    pi
  CHARACTER (100)  :: cbuff
!!!$.....................................................................
  pi = dacos(-1d0)

  lw_ref = .FALSE.
  lam_ref = 1d0 

!!!$     'datei' oeffnen
  fetxt = datei

  errnr = 1
  OPEN(kanal,file=TRIM(fetxt),status='old',err=999)

  errnr = 3
  IF (lnsepri) THEN
     PRINT*,''
     CALL get_unit(ifp1)
     PRINT*,'iseedpri / std. err.',iseedpri,modl_stdn
     OPEN (ifp1,FILE='tmp.mynoise_mprior_rho',STATUS='replace')
     WRITE (*,'(A)',advance='no')' Initializing noise '
     ALLOCATE (rnd_r(elanz))
     CALL Random_Init(iseedpri)
     DO i=1,elanz
        rnd_r(i) = Random_Gauss()
     END DO
     IF (.NOT.ldc) THEN
        CALL get_unit(ifp2)
        OPEN (ifp2,FILE='tmp.mynoise_mprior_phase',STATUS='replace')
        ALLOCATE (rnd_p(elanz))
        CALL Random_Init(-iseedpri)
        DO i=1,elanz
           rnd_p(i) = Random_Gauss()
        END DO
     END IF
     PRINT*,''
  END IF

!!!$     Anzahl der Messwerte lesen
!!!$ also check if we may use individual errors or not
!!! Failsafe read in one line

  READ(kanal,*,END=11,err=11) idum,lw_ref,lam_ref
  PRINT*,'Successful reference options'
  GOTO 12
11 lam_ref = 1d0 
  BACKSPACE (kanal)

12 IF (lw_ref) THEN

     PRINT*,'---> Reference model regularization '
     IF (lam_ref <= EPSILON(REAL(lam_ref))) THEN
        lam_ref = 0d0 ! so, what again is zero ??
     END IF
     PRINT*,'     lambda ref (factor)=',lam_ref
  END IF

!!!$     Ggf. Fehlermeldung
  IF (idum.NE.elanz) THEN
     fetxt = ' '
     errnr = 47
     GOTO 1000
  END IF
  


!!!$     Betrag und Phase (in mrad) des komplexen Widerstandes einlesen
  DO i=1,elanz
     
     IF (lw_ref) THEN
!!$        IF (ldc) THEN
!!$           READ(kanal,*,END=1001,err=1002) bet,pha,bet_ref,w_ref_re(mnr(i))
!!$           pha_ref = 0d0
!!$        ELSE
        READ(kanal,*,END=1001,err=1002) bet,pha,bet_ref,eps_r,pha_ref,eps_p
        w_ref_re(mnr(i)) = eps_r * eps_r
        w_ref_im(mnr(i)) = eps_p * eps_p
        pha_ref = 1d-3 * pha_ref ! mRad -> Rad
!!$        END IF
     ELSE
        READ(kanal,*,END=1001,err=1000) bet,pha
     END IF

     pha      = 1d-3*pha

     IF (lprior) THEN 
        !     set prior model ...             
        IF (bet0 <= 0. .OR. &
             (.NOT.ldc.AND.dabs(pha0).GT.1d3*pi)) THEN
           fetxt = 'starting model incorrect '
           errnr = 91
           GOTO 999
        END IF

        IF (bet > 0d0) THEN  
           !     TODO: meaningful phase check.. 
           IF (lnsepri) THEN
              eps_r = 1d-2*modl_stdn * bet
              WRITE(ifp1,'(3(G14.4,1X))',ADVANCE='no')rnd_r(i),eps_r,bet
              bet = bet + rnd_r(i) * eps_r
              WRITE(ifp1,'(G14.4)')bet
              IF (.NOT. ldc) THEN
                 eps_p = 1d-2*modl_stdn*dabs(pha)
                 WRITE(ifp2,'(3(G14.4,1X))',ADVANCE='no')rnd_p(i),eps_p,pha
                 pha = pha + rnd_p(i) * eps_p
                 WRITE(ifp2,'(G14.4)')pha
              END IF
           END IF
           sigma(i) = dcmplx(dcos(pha)/bet,-dsin(pha)/bet)
!!!$    hier koennte mittelung passieren
           m0(mnr(i)) = cdlog(sigma(i))
        ELSE                
           !     or let it stay at zero and take background cond
           sigma(i) = dcmplx( dcos(pha0/1d3)/bet0 , -dsin(pha0/1d3)/bet0 )
           m0(mnr(i)) = DCMPLX(0d0)
        END IF

     ELSE

!!!$     Ggf. Fehlermeldung
        IF (bet.LT.1d-12) THEN
           fetxt = ' '
           errnr = 11
           GOTO 999
        END IF

!!!$     Komplexe Leitfaehigkeit bestimmen
        sigma(i) = dcmplx(dcos(pha)/bet,-dsin(pha)/bet)

     END IF

!!!$ >> RM ref model regu
     IF (.NOT.lw_ref) CYCLE ! next i 
     IF (w_ref_re(mnr(i)) > EPSILON(eps_r)) THEN
!!!$ assign m_ref = \ln(|\sigma|) + i \phi(\sigma) 
!!!$              = -\ln(|\rho|) - i\phi(\rho) 
!!!$              = - (\ln(|\rho|)+i\phi(\rho)
!!!$ for \phi(z) = Im(z) / Re(z), z\inC

        m_ref(mnr(i)) = -DCMPLX(DLOG(bet_ref),pha_ref)
        !        PRINT*,'assign',mnr(i),wref(mnr(i))
     ELSE
        m_ref(mnr(i)) = DCMPLX(0d0)
!        PRINT*,'dont assign',mnr(i)
     END IF
!!!$ << RM ref model regu

  END DO

!!!$     'datei' schliessen
  CLOSE(kanal)

  IF (lnsepri) THEN
     CLOSE (ifp1)
     IF (.NOT.ldc) CLOSE(ifp2)
  END IF

  errnr = 0
  RETURN

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

999 RETURN

1000 CLOSE(kanal)
  RETURN

1001 CLOSE(kanal)
  errnr = 2
  RETURN

1002 CLOSE(kanal)
  fetxt = 'reference assignment failed'
  errnr = 2
  RETURN

END SUBROUTINE rsigma

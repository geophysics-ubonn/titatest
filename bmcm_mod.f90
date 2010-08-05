MODULE bmcm_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!$ Collection of subroutines to calculate linearized model 
!!!$ uncertainty and resolution matrices for ERT (DC) and EIT (IP)
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
  USE alloci , ONLY : sens,sensdc,smatm,nachbar,&
       ata_dc,ata_reg_dc,cov_m_dc,ata,ata_reg,cov_m
  USE femmod , ONLY : ldc
  USE elemmod, ONLY : smaxs
  USE invmod , ONLY : lip,wmatd,wdfak
  USE errmod , ONLY : errnr,fetxt
  USE konvmod , ONLY : ltri,lgauss,lam,nx,nz,mswitch,lcov2,lres
  USE modelmod , ONLY : manz
  USE datmod , ONLY : nanz
  USE errmod, ONLY : errnr,fetxt,fprun
  USE sigmamod , ONLY : sigma
  USE pathmod

  IMPLICIT none

  PUBLIC :: buncert
!!!$ this sub controls uncertainty caculation

  PRIVATE :: bata_dc
!!!$ A^TC_d^-1A  -> ata_dc
  PRIVATE :: bata_reg_dc 
!!!$ ata_dc + lam*C_m^-1 -> ata_reg_dc
  PRIVATE :: bmcm_dc ! inversion is done by cholesky or gauss
!!!$ ata_reg_dc^-1   -> cov_m_dc
  PRIVATE :: bres_dc ! solution with MATMUL
!!!$ cov_m_dc * ata_dc -> ata_reg_dc 
  PRIVATE :: bmcm2_dc ! using MATMUL
!!!$ ata_reg_dc * cov_m_dc -> ata_dc

  PRIVATE :: bata
!!!$ A^TC_d^-1A  -> ata
  PRIVATE :: bata_reg
!!!$ ata + lam*C_m^-1 -> ata_reg
  PRIVATE :: bmcm ! inversion
!!!$ ata_reg^-1   -> cov_m
  PRIVATE :: bres ! MATMUL
!!!$ cov_m * ata -> ata_reg 
  PRIVATE :: bmcm2 ! MATMUL
!!!$ ata_reg * cov_m -> ata

CONTAINS

  SUBROUTINE buncert (kanal,lamalt)

!!!$
!!!$ This sub is the control unit of the smatm calculation
!!!$
    INTEGER (KIND = 4 ),INTENT(IN) :: kanal ! kanal is the io unit number
    REAL (KIND(0D0)),INTENT(IN)    :: lamalt ! lambda of the last iteration
!!! for tic_toc
    INTEGER (KIND = 4 )            :: c1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!$    get time
    CALL TIC(c1)

    WRITE(*,'(a)')'Calculating model uncertainty..'
    WRITE (fprun,'(a)')'Calculating model uncertainty..'
    lam = lamalt
    WRITE (*,'(/a,G10.3,a/)')'take current lambda ?',lam,&
         ACHAR(9)//':'//ACHAR(9)
    IF (BTEST(mswitch,5)) THEN 
       READ (*,*)fetxt
       IF (fetxt/='')READ(fetxt,*)lam
       WRITE (*,*)'Set lambda to ',lam
    END IF

    WRITE (fprun,*)'Taking lam=',lam

    WRITE(*,'(a)')ACHAR(13)//&
         'calculating MCM_1 = (A^TC_d^-1A + C_m^-1)^-1'
    WRITE(fprun,'(a)')'MCM_1 = (A^TC_d^-1A + C_m^-1)^-1'

    IF (ldc) THEN

       ALLOCATE (ata_dc(manz,manz),STAT=errnr)
       IF (errnr /= 0) THEN
          errnr = 97
          RETURN
       END IF
       ata_dc = 0D0

       fetxt = ramd(1:lnramd)//slash(1:1)//'ata.diag'
       CALL bata_dc(kanal) ! A^TC_d^-1A -> ata_dc
       if (errnr.ne.0) RETURN

       ALLOCATE (ata_reg_dc(manz,manz),STAT=errnr)            
       IF (errnr /= 0) THEN
          errnr = 97
          RETURN
       END IF

       fetxt = ramd(1:lnramd)//slash(1:1)//'ata_reg.diag'
       CALL bata_reg_dc(kanal) !ata_dc + C_m^-1(lam) -> ata_reg_dc
       if (errnr.ne.0) RETURN

       ALLOCATE (cov_m_dc(manz,manz))
       IF (errnr/=0) THEN
          WRITE (*,'(/a/)')'Allocation problem MCM_1 in bmcmdc'
          errnr = 97
          RETURN
       END IF

       fetxt = ramd(1:lnramd)//slash(1:1)//'cov1_m.diag'
       CALL bmcm_dc(kanal,lgauss) ! (A^TC_d^-1A + C_m^-1)^-1   -> cov_m_dc
       if (errnr.ne.0) RETURN

    ELSE

       ALLOCATE (ata(manz,manz),STAT=errnr)
       IF (errnr /= 0) THEN
          errnr = 97
          RETURN
       END IF
       ata = DCMPLX(0.)

       fetxt = ramd(1:lnramd)//slash(1:1)//'ata.diag'
       CALL bata(kanal)    ! A^TC_d^-1A -> ata
       if (errnr.ne.0) RETURN

       ALLOCATE (ata_reg(manz,manz),STAT=errnr)
       IF (errnr /= 0) THEN
          errnr = 97
          RETURN
       END IF

       fetxt = ramd(1:lnramd)//slash(1:1)//'ata_reg.diag'
       CALL bata_reg(kanal) !A^TC_d^-1A + C_m^-1(lam) -> ata_reg
       if (errnr.ne.0) RETURN

       ALLOCATE (cov_m(manz,manz))
       IF (errnr/=0) THEN
          WRITE (*,'(/a/)')'Allocation problem MCM_1 in bmcmdc'
          errnr = 97
          RETURN
       END IF

       fetxt = ramd(1:lnramd)//slash(1:1)//'cov1_m.diag'
       CALL bmcm(kanal,lgauss)    ! (A^TC_d^-1A + C_m^-1)^-1   -> cov_m
       if (errnr.ne.0) RETURN

    END IF

    IF (lres) THEN

       fetxt = ramd(1:lnramd)//slash(1:1)//'res_m.diag'

       WRITE(*,'(a)')ACHAR(13)//'calculating RES = '//&
            '(A^TC_d^-1A + C_m^-1)^-1A^TC_d^-1A'
       WRITE(fprun,'(a)')'RES = (A^TC_d^-1A + C_m^-1)^-1A^TC_d^-1A'

       IF (ldc) THEN
          CALL bres_dc(kanal)
          if (errnr.ne.0) RETURN
       ELSE
          CALL bres(kanal)
          if (errnr.ne.0) RETURN
       END IF

       IF (lcov2) THEN

          fetxt = ramd(1:lnramd)//slash(1:1)//'cov2_m.diag'

          WRITE(*,'(a)')ACHAR(13)//'calculating MCM_2 = '//&
               '(A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A (A^TC_d^-1A + C_m^-1)^-1'
          WRITE(fprun,'(a)')'MCM_2 = (A^TC_d^-1A + C_m^-1)'//&
               '^-1 A^TC_d^-1A (A^TC_d^-1A + C_m^-1)^-1'

          IF (ldc) THEN
             CALL bmcm2_dc(kanal)
             if (errnr.ne.0) RETURN
          ELSE
             CALL bmcm2(kanal)
             if (errnr.ne.0) RETURN
          END IF

       END IF              ! lcov2
    END IF                 ! lres

!!! $ may be we ant to write out not only main diagonals before freeing

    IF (ALLOCATED (ata)) DEALLOCATE (ata)
    IF (ALLOCATED (ata_dc)) DEALLOCATE (ata_dc)
    IF (ALLOCATED (ata_reg)) DEALLOCATE (ata_reg)
    IF (ALLOCATED (ata_reg_dc)) DEALLOCATE (ata_reg_dc)
    IF (ALLOCATED (cov_m)) DEALLOCATE (cov_m)
    IF (ALLOCATED (cov_m_dc)) DEALLOCATE (cov_m_dc)
    
    fetxt = 'uncertainty calculation time'
    CALL TOC(c1,fetxt)

  END SUBROUTINE buncert

  SUBROUTINE bata_dc(kanal)
!!!$    
!!!$    Unterprogramm berechnet ATC_d^-1A
!!!$
!!!$    Copyright Andreas Kemna 2009
!!!$    Created by  Roland Martin                            02-Nov-2009
!!!$    
!!!$    Last changed      RM                                   20-Feb-2010
!!!$
!!!$.....................................................................
!!!$   PROGRAMMINTERNE PARAMETER:
!!!$   Hilfsvariablen 
    INTEGER                       :: kanal
    INTEGER                       :: i,j,k
    REAL,DIMENSION(:),ALLOCATABLE :: dig
    REAL                          :: dig_min,dig_max
!!!$....................................................................

!!!$  A^TC_d^-1A

    errnr = 1
    OPEN (kanal,file=TRIM(fetxt),status='replace',err=999)
    errnr = 4

    ALLOCATE(dig(manz))

    ata_dc= 0.
    DO k=1,manz
       write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//&
            'ATC_d^-1A/ ',REAL( k * (100./manz)),'%'
       DO j=k,manz ! fills upper triangle (k,j)
          DO i=1,nanz
             ata_dc(k,j) = ata_dc(k,j) + sensdc(i,k) * & 
                  sensdc(i,j) * wmatd(i) * DBLE(wdfak(i))
          END DO
          ata_dc(j,k) = ata_dc(k,j) ! fills lower triangle (k,j)
       END DO
       dig(k) = REAL(ata_dc(k,k))
    END DO
!!!$    get min/max
    dig_min = MINVAL(dig)
    dig_max = MAXVAL(dig)
!!!$    write out 
    WRITE (kanal,*)manz
    DO i=1,manz
       WRITE (kanal,*)LOG10(dig(i)/dig_max),dig(i)
    END DO
    WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
    WRITE (*,*)'Max/Min:',dig_max,'/',dig_min
    CLOSE(kanal)

    DEALLOCATE (dig)

    errnr = 0
999 RETURN

  END SUBROUTINE bata_dc

  SUBROUTINE bata_reg_dc(kanal)
!!!$    
!!!$    Unterprogramm berechnet ATC_d^-1A+lam*C_m
!!!$    Fuer beliebige Triangulierung
!!!$    
!!!$    Copyright Andreas Kemna
!!!$    Erstellt von Roland Martin                           02-Nov-2009
!!!$    
!!!$    Last changed      RM                               31-Mar-2010
!!!$    
!!!$.....................................................................
!!!$   PROGRAMMINTERNE PARAMETER:
    INTEGER                        :: kanal ! io number
!!!$   Hilfsvariablen 
    INTEGER                        :: i,j,k
    REAL,DIMENSION(:),ALLOCATABLE  :: dig
    REAL                           :: dig_min,dig_max
!!!$.....................................................................

!!!$  A^TC_d^-1A+lamC_m

    errnr = 1
    open(kanal,file=TRIM(fetxt),status='replace',err=999)
    errnr = 4

    ALLOCATE(dig(manz))

    IF (ltri == 0) THEN
       DO j=1,manz
          write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)// &
               'ATC_d^-1A+reg/ ',REAL( j * (100./manz)),'%'
          DO i=1,j            ! lower triangle

             IF (i == j) THEN
                ata_reg_dc(i,j) = ata_dc(i,j) + lam * smatm(i,1)

                IF (i+1 < manz) ata_reg_dc(i+1,j) = ata_dc(i+1,j) + &
                     lam * smatm(i+1,2) ! nebendiagonale in x richtung
                IF (i+nx < manz) ata_reg_dc(i+nx,j) = ata_dc(i+nx,j) + &
                     lam * smatm(i+nx,3) ! nebendiagonale in z richtung
             ELSE
                ata_reg_dc(i,j) = ata_dc(i,j) ! only aTa
             END IF

             ata_reg_dc(j,i) = ata_reg_dc(i,j) ! upper triangle 

          END DO

          dig(j) = REAL(ata_reg_dc(j,j))
       END DO

    ELSE IF (ltri == 3.OR.ltri == 4) THEN

       ata_reg_dc= ata_dc

       DO i=1,manz
          dig (i) = ata_dc(i,i) + lam * smatm(i,1)
          ata_reg_dc(i,i) = dig(i)
       END DO

    ELSE IF (ltri == 1.OR.ltri == 2.OR.(ltri > 4 .AND. ltri < 15)) THEN

       DO j=1,manz

          write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)// &
               'ATC_d^-1A+reg/ ',REAL( j * (100./manz)),'%'

          DO i=1,manz

             IF (i==j) THEN   ! sparse C_m

                ata_reg_dc(i,j) = ata_dc(i,j) + lam * smatm(i,smaxs+1)

                DO k=1,nachbar(i,smaxs+1)

                   IF (nachbar(i,k) /= 0) ata_reg_dc(nachbar(i,k),j) = &
                        ata_dc(nachbar(i,k),j) + lam * smatm(i,k)

                END DO

             ELSE

                ata_reg_dc(i,j) = ata_dc(i,j)

             END IF

          END DO

          dig(j) = REAL(ata_reg_dc(j,j))

       END DO

    ELSE IF (ltri == 15) THEN

       ata_reg_dc= ata_dc+ smatm ! for full C_m..w

       DO j=1,manz
          dig(j) = REAL(ata_reg_dc(j,j))
       END DO

    END IF

    dig_min = MINVAL(dig)
    dig_max = MAXVAL(dig)

    WRITE (kanal,*)manz
    DO i=1,manz
       WRITE (kanal,*)LOG10(dig(i)),dig(i)
    END DO

    WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
    WRITE (*,*)'Max/Min:',dig_max,'/',dig_min

    CLOSE(kanal)

    DEALLOCATE (dig)

    errnr = 0
999 RETURN

  END SUBROUTINE bata_reg_dc

  SUBROUTINE bmcm_dc(kanal,ols)
!!!$    
!!!$    Unterprogramm berechnet (einfache) Modell Kovarianz Matrix
!!!$    (A^TC_d^-1A + C_m^-1)^-1
!!!$    Fuer beliebige Triangulierung
!!!$    
!!!$    Copyright Andreas Kemna 2009
!!!$    Created by  Roland Martin                            02-Nov-2009
!!!$    
!!!$    Last changed      RM                                   20-Feb-2010
!!!$    
!!!$.....................................................................
!!!$   PROGRAMMINTERNE PARAMETER:
!!!$   Hilfsvariablen 
    INTEGER                                    :: i,kanal,j,c1
    REAL(KIND(0D0)),DIMENSION(:,:),ALLOCATABLE :: work
    REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: dig,dig2
    REAL(KIND(0D0))                            :: dig_min,dig_max
    LOGICAL,INTENT(IN),OPTIONAL                :: ols
    CHARACTER(80)                              :: csz
!!!$....................................................................

!!!$  invert (A^TC_d^-1A + C_m^-1)

    errnr = 1

    ALLOCATE (work(manz,manz),STAT=errnr)
    IF (errnr/=0) THEN
       fetxt = 'Allocation problem WORK in bmcm'
       errnr = 97
       RETURN
    END IF
    ALLOCATE (dig(manz),dig2(manz)) !prepare to write out main diagonal

!!!$    get time
    CALL TIC(c1)

    cov_m_dc= ata_reg_dc

    IF (.NOT.PRESENT(ols).OR..NOT.ols) THEN !default

       WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//'Factorization...'

       CALL CHOLD(cov_m_dc,dig,manz,errnr)
       IF (errnr /= 0) THEN
          fetxt='CHOLD mcm :: matrix not pos definite..'
          PRINT*,'Zeile(',abs(errnr),')'
          errnr = 108
          RETURN
       END IF
       WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//'Inverting...'

       CALL LINVD(cov_m_dc,dig,manz)

       WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//'Filling lower Cov...'

       DO i= 1 , manz

          WRITE (*,'(A,1X,F6.2,A)',ADVANCE='no')ACHAR(13)//ACHAR(9)&
               //ACHAR(9)//ACHAR(9)//'/ ',REAL( i * (100./manz)),'%'

          DO j = 1 , i - 1

             cov_m_dc(i,j) = cov_m_dc(j,i)

          END DO
       END DO

    ELSE
       WRITE (*,'(a)',ADVANCE='no')'Inverting Matrix (Gauss elemination)'

       CALL gauss_dble(cov_m_dc,manz,errnr)

       IF (errnr /= 0) THEN
          fetxt = 'error matrix inverse not found'
          PRINT*,'Zeile::',cov_m_dc(abs(errnr),:)
          PRINT*,'Spalte::',cov_m_dc(:,abs(errnr))
          errnr = 108
          RETURN
       END IF
    END IF

    csz = 'solution time'
    CALL TOC(c1,csz)

    work = MATMUL(cov_m_dc,ata_reg_dc)

    DO i=1,manz
       dig(i) = cov_m_dc(i,i)
       dig2(i) = work(i,i)
    END DO

    dig_min = MINVAL(dig)
    dig_max = MAXVAL(dig)

    errnr = 1
    OPEN (kanal,file=TRIM(fetxt)//'_re',status='replace',err=999)
    errnr = 4

    WRITE (kanal,*)manz
    DO i=1,manz
       WRITE (kanal,*)LOG10(SQRT(dig(i))),dig2(i)
    END DO

    WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
    WRITE (*,*)'Max/Min:',dig_max,'/',dig_min

    CLOSE(kanal)

    DEALLOCATE (dig,dig2,work)

    errnr = 0
999 RETURN

  END SUBROUTINE bmcm_dc

  SUBROUTINE bres_dc(kanal)
!!!$    
!!!$    Unterprogramm berechnet Aufloesungsmatrix
!!!$    Fuer beliebige Triangulierung
!!!$    RES = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A
!!!$    wobei (A^TC_d^-1A + C_m^-1) bereits invertiert wurde (cov_m_dc)
!!!$
!!!$    Copyright Andreas Kemna 2009
!!!$
!!!$    Created by Roland Martin                            02-Nov-2009
!!!$    
!!!$    Last changed      RM                             20-Feb-2010
!!!$.....................................................................
!!!$   PROGRAMMINTERNE PARAMETER:
!!!$   Hilfsvariablen 
    INTEGER                                      :: i,kanal
    REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE     :: dig
    REAL(KIND(0D0))                              :: dig_min,dig_max
!!!$.....................................................................

!!!$  cal!!!$RES = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A

    errnr = 1
    open(kanal,file=fetxt,status='replace',err=999)
    errnr = 4

    ata_reg_dc= MATMUL(cov_m_dc,ata_dc) ! that's it...

    ALLOCATE (dig(manz)) !prepare to write out main diagonal

    DO i=1,manz
       dig(i) = ata_reg_dc(i,i)
    END DO

    dig_min = MINVAL(dig)
    dig_max = MAXVAL(dig)

    WRITE (kanal,*)manz
    DO i=1,manz
       WRITE (kanal,*)ABS(dig(i)),LOG10(dig(i))
    END DO

    WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
    WRITE (*,*)'Max/Min:',dig_max,'/',dig_min

    CLOSE(kanal)

    DEALLOCATE (dig)

    errnr = 0
999 RETURN

  END SUBROUTINE bres_dc

  SUBROUTINE bmcm2_dc(kanal)
!!!$    
!!!$    Unterprogramm berechnet Modellkovarianz nach Fehlerfortpflanzung
!!!$    MCM = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A (A^TC_d^-1A + C_m^-1)^-1
!!!$    Fuer beliebige Triangulierung
!!!$!!!$    Copyright Andreas Kemna/Roland Martin 2009
!!!$    erstellt von Roland Martin                               02-Nov-2009
!!!$    
!!!$    Last changed      RM                                   23-Nov-2009
!!!$    
    !!!$.........................................................................
    USE alloci
    USE datmod
    USE invmod
    USE modelmod
    USE errmod

    IMPLICIT none

!!!$.....................................................................
!!!$   PROGRAMMINTERNE PARAMETER:
!!!$   Hilfsvariablen 
    INTEGER                                      :: kanal
    INTEGER                                      :: i,j,k
    REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE     :: dig
    REAL(KIND(0D0))                              :: dig_min,dig_max
!!!$.....................................................................
!!!$   vorher wurde schon 
!!!$   ata_reg_dc= (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A berechnet
!!!$   cov_m_dc= (A^TC_d^-1A + C_m^-1)^-1


    errnr = 1
    open(kanal,file=fetxt,status='replace',err=999)
    errnr = 4

    ata_dc = MATMUL (ata_reg_dc,cov_m_dc)

    ALLOCATE (dig(manz))

    DO i=1,manz
       dig(i) = ata_dc(i,i)
    END DO

    dig_min = MINVAL(dig) 
    dig_max = MAXVAL(dig)

    WRITE (kanal,*)manz
    DO i=1,manz
       WRITE (kanal,*)LOG10(SQRT(ABS(dig(i)))),dig(i)
    END DO

    WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
    WRITE (*,*)'Max/Min:',dig_max,'/',dig_min

    CLOSE(kanal)

    DEALLOCATE (dig)

    errnr = 0

999 RETURN

  END SUBROUTINE bmcm2_dc
  SUBROUTINE bata(kanal)
!!!$    
!!!$    Unterprogramm berechnet ATA=A^T C_d^{-1} A (complex)
!!!$    
!!!$    Copyright Andreas Kemna 2009
!!!$    Created by  Roland Martin                            02-Nov-2009
!!!$    
!!!$    Last changed      RM                                   20-Feb-2010
!!!$    !!!$.........................................................................

    USE alloci
    USE datmod
    USE invmod
    USE modelmod
    USE errmod

    IMPLICIT none

!!!$.....................................................................
!!!$   PROGRAMMINTERNE PARAMETER:
!!!$   Hilfsvariablen 
    INTEGER                       :: kanal
    INTEGER                       :: i,j,k
    REAL,DIMENSION(:),ALLOCATABLE :: dig
    REAL                          :: dig_min,dig_max
!!!$.....................................................................

!!!$  A^TC_d^-1A

    errnr = 1
    OPEN(kanal,file=TRIM(fetxt),status='replace',err=999)
    errnr = 4

    ALLOCATE(dig(manz))

    DO k=1,manz
       write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//&
            'ATC_d^-1A/ ',REAL( k * (100./manz)),'%'
       DO j=k,manz ! fills upper triangle (k,j) 
          DO i=1,nanz
             ata(k,j) = ata(k,j) + DCONJG(sens(i,k)) * &
                  sens(i,j) * wmatd(i) * DBLE(wdfak(i))
          END DO
          ata(j,k) = ata(k,j) ! fills lower triangle (j,k) 
       END DO
       dig(k) = REAL(ata(k,k))
    END DO

!!!$    write out 
    dig_min = MINVAL(dig)
    dig_max = MAXVAL(dig)
    WRITE (kanal,*)manz
    DO i=1,manz
       WRITE (kanal,*)LOG10(dig(i)/dig_max),dig(i)
    END DO
    WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
    WRITE (*,*)'Max/Min:',dig_max,'/',dig_min
    CLOSE(kanal)

    DEALLOCATE (dig)

    errnr = 0
999 RETURN

  END SUBROUTINE bata

  SUBROUTINE bata_reg(kanal)
!!!$    
!!!$    Unterprogramm berechnet ATC_d^-1A+lam*C_m
!!!$    Fuer beliebige Triangulierung
!!!$    
!!!$    Copyright Andreas Kemna
!!!$    Created by  Roland Martin                            02-Nov-2009
!!!$    
!!!$    Last changed      RM                                   20-Feb-2010
!!!$    
    !!!$.........................................................................

    USE alloci
    USE invmod
    USE modelmod
    USE elemmod
    USE errmod
    USE konvmod

    IMPLICIT none

!!!$.....................................................................
!!!$   PROGRAMMINTERNE PARAMETER:
    INTEGER                       :: kanal ! io number
!!!$   Hilfsvariablen 
    INTEGER                       :: i,j,k
    REAL,DIMENSION(:),ALLOCATABLE :: dig
    REAL                          :: dig_min,dig_max
!!!$.....................................................................

!!!$  A^TC_d^-1A+lamC_m

    errnr = 1
    OPEN(kanal,file=TRIM(fetxt),status='replace',err=999)
    errnr = 4

    ALLOCATE(dig(manz))

    IF (ltri == 0) THEN
       DO j=1,manz
          write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//&
               'ATC_d^-1A+reg/ ',REAL( j * (100./manz)),'%'
          DO i=1,j ! lower triangle

             IF (i == j) THEN
                ata_reg(i,j) = ata(i,j) + DCMPLX(lam * smatm(i,1))

                IF (i+1 < manz) ata_reg(i+1,j) = ata(i+1,j) + &
                     DCMPLX(lam * smatm(i+1,2)) ! nebendiagonale in x richtung
                IF (i+nx < manz) ata_reg(i+nx,j) = ata(i+nx,j) + &
                     DCMPLX(lam * smatm(i+nx,3)) ! nebendiagonale in z richtung
             ELSE
                ata_reg(i,j) = ata(i,j) ! only aTa
             END IF

             ata_reg(j,i) = ata_reg(i,j) ! upper triangle 

          END DO
          dig(j) = REAL(ata_reg(j,j))
       END DO

    ELSE IF (ltri == 3.OR.ltri == 4) THEN
       ata_reg = ata
       DO i=1,manz
          ata_reg(i,i) = ata(i,j) + DCMPLX(lam * smatm(i,1))
          dig(i) = REAL(ata_reg(i,i))
       END DO

    ELSE IF (ltri == 1.OR.ltri == 2.OR.(ltri > 4 .AND. ltri < 15)) THEN

       DO j=1,manz

          write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//&
               'ATC_d^-1A+reg/ ',REAL( j * (100./manz)),'%'

          DO i=1,manz

             IF (i==j) THEN   ! sparse C_m

                ata_reg(i,j) = ata(i,j) + lam * smatm(i,smaxs+1)

                DO k=1,nachbar(i,smaxs+1)
                   IF (nachbar(i,k) /= 0) ata_reg(nachbar(i,k),j) = &
                        ata(nachbar(i,k),j) + DCMPLX(lam * smatm(i,k))
                END DO

             ELSE

                ata_reg(i,j) = ata(i,j)

             END IF

          END DO

          dig(j) = REAL(ata_reg(j,j))

       END DO

    ELSE IF (ltri == 15) THEN

       ata_reg = ata + DCMPLX(lam) * DCMPLX(smatm) ! for full C_m..

       DO j=1,manz
          dig(j) = REAL(ata_reg(j,j))
       END DO

    END IF

!!!$    write out real and imaginary part
    dig_min = MINVAL(dig)
    dig_max = MAXVAL(dig)
    WRITE (kanal,*)manz
    DO i=1,manz
       WRITE (kanal,*)LOG10(dig(i)),dig(i)
    END DO
    WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
    WRITE (*,*)'Max/Min:',dig_max,'/',dig_min
    CLOSE(kanal)

    DEALLOCATE (dig)

    errnr = 0
999 RETURN

  END SUBROUTINE bata_reg

  SUBROUTINE bmcm(kanal,ols)
!!!$    
!!!$    Unterprogramm berechnet (einfache) Modell Kovarianz Matrix
!!!$    MCM = (A^TC_d^-1A + C_m^-1)^-1
!!!$    Fuer beliebige Triangulierung
!!!$    
!!!$    Copyright Andreas Kemna
!!!$    Created by  Roland Martin                      02-Nov-2009
!!!$    
!!!$    Last changed      RM                             30-Mar-2010
!!!$    
!!!$....................................................................
!!!$   PROGRAMMINTERNE PARAMETER:
!!!$   Hilfsvariablen 
    INTEGER                                     :: i,kanal,j,c1
    COMPLEX(KIND(0D0)),DIMENSION(:,:),ALLOCATABLE :: work
    COMPLEX(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: dig,dig2
    REAL(KIND(0D0))                             :: dig_min,dig_max,p
    LOGICAL,INTENT(IN),OPTIONAL                 :: ols 
    CHARACTER(80)                               :: csz
!!!$....................................................................
!!!$  invert A^TC_d^-1A + C_m^-1
    errnr = 1

    ALLOCATE (work(manz,manz),STAT=errnr)
    IF (errnr/=0) THEN
       WRITE (*,'(/a/)')'Allocation problem WORK in bmcm'
       errnr = 97
       RETURN
    END IF
    ALLOCATE (dig(manz),dig2(manz)) ! help 

!!!$    get time
    CALL TIC(c1)

    cov_m = ata_reg

    IF (.NOT.PRESENT(ols).OR..NOT.ols) THEN
       WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//'Factorization...'

       CALL CHOLZ(cov_m,dig,manz,errnr)
       IF (errnr /= 0) THEN
          PRINT*,'Zeile::',cov_m(abs(errnr),:)
          PRINT*,'Spalte::',cov_m(:,abs(errnr))
          errnr = 108
          RETURN
       END IF
       WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//'Inverting...'
       CALL LINVZ(cov_m,dig,manz)

       WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//'Filling lower Cov...'

       DO i= 1 , manz

          WRITE (*,'(A,1X,F6.2,A)',ADVANCE='no')ACHAR(13)//ACHAR(9)&
               //ACHAR(9)//ACHAR(9)//'/ ',REAL( i * (100./manz)),'%'

          DO j = 1 , i - 1

             cov_m(i,j) = cov_m(j,i)

          END DO
       END DO
    ELSE
       WRITE (*,'(a)',ADVANCE='no')'Inverting Matrix (Gauss elemination)'
       CALL gauss_cmplx(cov_m,manz,errnr)

       IF (errnr /= 0) THEN
          PRINT*,'Zeile::',cov_m(abs(errnr),:)
          PRINT*,'Spalte::',cov_m(:,abs(errnr))
          errnr = 108
          RETURN
       END IF

    END IF

    csz = 'solution time'
    CALL TOC(c1,csz)

    work = MATMUL(cov_m,ata_reg)

    DO i=1,manz
       dig(i) = cov_m(i,i)
       dig2(i) = work(i,i)
    END DO

!!!$    write out real and imaginary part
    errnr = 1
    OPEN(kanal,file=TRIM(fetxt)//'_re',status='replace',err=999)
    errnr = 4
    dig_min = MINVAL(DBLE(dig))
    dig_max = MAXVAL(DBLE(dig))
    WRITE (kanal,*)manz
    DO i=1,manz
       WRITE (kanal,*)LOG10(SQRT(DBLE(dig(i)))),DBLE(dig2(i))
    END DO
    WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
    WRITE (*,*)'Max/Min(Re):',dig_max,'/',dig_min
    CLOSE(kanal)

    errnr = 1
    OPEN(kanal,file=TRIM(fetxt)//'_im',status='replace',err=999)
    errnr = 4
    dig_min = MINVAL(DIMAG(dig))
    dig_max = MAXVAL(DIMAG(dig))
    WRITE (kanal,*)manz
    DO i=1,manz
       p = DBLE(1d3*DATAN2(DIMAG(dig(i)),DBLE(dig(i))))
       WRITE (kanal,*)DIMAG(dig(i)),p
    END DO
    WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
    WRITE (*,*)'Max/Min(Im):',dig_max,'/',dig_min
    CLOSE(kanal)


    DEALLOCATE (dig)
    IF (ALLOCATED(work)) DEALLOCATE (dig2,work)

    errnr = 0
999 RETURN
  END SUBROUTINE bmcm

  SUBROUTINE bres(kanal)
!!!$    
!!!$    Unterprogramm berechnet Aufloesungsmatrix
!!!$    Fuer beliebige Triangulierung und Complex
!!!$    RES = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A
!!!$    wobei (A^TC_d^-1A + C_m^-1) bereits invertiert wurde (cov_m)
!!!$    
!!!$    Copyright Andreas Kemna 2009
!!!$    Created by  Roland Martin                      02-Nov-2009
!!!$    
!!!$    Last changed      RM                             30-Mar-2010
!!!$    
!!!$....................................................................
!!!$   PROGRAMMINTERNE PARAMETER:
!!!$   Hilfsvariablen 
    INTEGER                                     :: i,kanal
    COMPLEX(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: dig
    REAL(KIND(0D0))                             :: dig_min,dig_max
!!!$   switch solv ordinary linear system or not^^
!!!$....................................................................

!!!$  calculate  RES = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A

    ata_reg = MATMUL(cov_m,ata)

    ALLOCATE (dig(manz))

    DO i=1,manz
       dig(i) = ata_reg(i,i)
    END DO

!!!$    write out real and imaginary part
    errnr = 1
    OPEN(kanal,file=TRIM(fetxt)//'_re',status='replace',err=999)
    errnr = 4
    dig_min = MINVAL(DBLE(dig))
    dig_max = MAXVAL(DBLE(dig))
    WRITE (kanal,*)manz
    DO i=1,manz
       WRITE (kanal,*)LOG10(DBLE(dig(i))),DBLE(dig(i))
    END DO
    WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
    WRITE (*,*)'Max/Min(Re):',dig_max,'/',dig_min
    CLOSE(kanal)

    errnr = 1
    OPEN(kanal,file=TRIM(fetxt)//'_im',status='replace',err=999)
    errnr = 4
    dig_min = MINVAL(DIMAG(dig))
    dig_max = MAXVAL(DIMAG(dig))
    WRITE (kanal,*)manz
    DO i=1,manz
       WRITE (kanal,*)LOG10(DIMAG(dig(i))),DIMAG(dig(i))
    END DO
    WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
    WRITE (*,*)'Max/Min(Im):',dig_max,'/',dig_min
    CLOSE(kanal)

    DEALLOCATE (dig)

    errnr = 0

999 RETURN

  END SUBROUTINE bres

  SUBROUTINE bmcm2(kanal)
!!!$    
!!!$    Unterprogramm berechnet Modellkovarianz nach 
!!!$    Fehlerfortpflanzung (Gubbins 2004)
!!!$    MCM = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A * 
!!!$          (A^TC_d^-1A + C_m^-1)^-1
!!!$    Fuer beliebige Triangulierung und Complex
!!!$
!!!$    Copyright Andreas Kemna 2009
!!!$
!!!$    Created by Roland Martin                      02-Nov-2009
!!!$    
!!!$    Last changed      RM                             30-Mar-2010
!!!$    
!!!$....................................................................
!!!$   PROGRAMMINTERNE PARAMETER:
!!!$   Hilfsvariablen 
    INTEGER                                     :: kanal
    INTEGER                                     :: i,j,k
    COMPLEX(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: dig
    COMPLEX(KIND(0D0))                          :: dum
    REAL(KIND(0D0))                             :: dig_min,dig_max,p
!!!$....................................................................
!!!$   vorher wurde schon 
!!!$   ata_reg = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A berechnet
!!!$   cov_m = (A^TC_d^-1A + C_m^-1)^-1

    ata = MATMUL (ata_reg,cov_m)

    ALLOCATE (dig(manz))

    DO i=1,manz
       dig(i) = ata(i,i)
    END DO

!!!$    write out real and imaginary part
    errnr = 1
    OPEN(kanal,file=TRIM(fetxt)//'_re',status='replace',err=999)
    errnr = 4
    dig_min = MINVAL(DBLE(dig))
    dig_max = MAXVAL(DBLE(dig))
    WRITE (kanal,*)manz
    DO i=1,manz
       WRITE (kanal,*)LOG10(SQRT(DBLE(dig(i)))),DBLE(dig(i))
    END DO
    WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
    WRITE (*,*)'Max/Min(Re):',dig_max,'/',dig_min
    CLOSE(kanal)

    errnr = 1
    OPEN(kanal,file=TRIM(fetxt)//'_im',status='replace',err=999)
    errnr = 4
    dig_min = MINVAL(DIMAG(dig))
    dig_max = MAXVAL(DIMAG(dig))
    WRITE (kanal,*)manz
    DO i=1,manz
       dum = dcmplx(1d0)/sigma(i)
       p = DBLE(1d3*datan2(dimag(dum),dble(dum)))
       WRITE (kanal,*)DIMAG(dig(i))*p,DIMAG(dig(i))
    END DO
    WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
    WRITE (*,*)'Max/Min(Im):',dig_max,'/',dig_min
    CLOSE(kanal)


    DEALLOCATE (dig)

    errnr = 0
999 RETURN

  END SUBROUTINE bmcm2

END MODULE bmcm_mod

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
       ata,ata_reg,cov_m
  USE femmod , ONLY : ldc
  USE elemmod, ONLY : smaxs,espx,espy
  USE invmod , ONLY : lip,wmatd,wdfak
  USE errmod , ONLY : errnr,fetxt
  USE konvmod , ONLY : ltri,lgauss,lam,nx,nz,mswitch,lcov2,lres,lverb,lverb_dat
  USE modelmod , ONLY : manz
  USE datmod , ONLY : nanz
  USE errmod, ONLY : errnr,fetxt,fprun
  USE sigmamod , ONLY : sigma
  USE pathmod
  USE ompmod

  IMPLICIT none

  PUBLIC :: buncert
!!!$ this sub controls uncertainty caculation

  PRIVATE :: bata
!!!$ A^HC_d^-1A  -> ata
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
    IF (BTEST(mswitch,6)) THEN 
       READ (*,*)fetxt
       IF (fetxt/='')READ(fetxt,*)lam
       WRITE (*,*)'No, set lambda to ',lam
    ELSE
       WRITE (*,*)' Yes'
    END IF

    WRITE (fprun,*)'Taking lam=',lam

    WRITE(*,'(a)')ACHAR(13)//&
         'calculating MCM_1 = (A^TC_d^-1A + C_m^-1)^-1'
    WRITE(fprun,'(a)')'MCM_1 = (A^TC_d^-1A + C_m^-1)^-1'

    ALLOCATE (ata(manz,manz),STAT=errnr)
    IF (errnr /= 0) THEN
       errnr = 97
       RETURN
    END IF
    ata = 0D0
    
    fetxt = ramd(1:lnramd)//slash(1:1)//'ata.diag'
    CALL bata(kanal) ! A^TC_d^-1A -> ata
    if (errnr.ne.0) RETURN
    
    ALLOCATE (ata_reg(manz,manz),STAT=errnr)
    IF (errnr /= 0) THEN
       errnr = 97
       RETURN
    END IF
    ata_reg = 0D0
    
    fetxt = ramd(1:lnramd)//slash(1:1)//'ata_reg.diag'
    CALL bata_reg(kanal) !ata_dc + C_m^-1(lam) -> ata_reg_dc
    if (errnr.ne.0) RETURN
    
    ALLOCATE (cov_m(manz,manz))
    IF (errnr/=0) THEN
       WRITE (*,'(/a/)')'Allocation problem MCM_1 in bmcmdc'
       errnr = 97
       RETURN
    END IF
    cov_m = 0D0
    
    fetxt = ramd(1:lnramd)//slash(1:1)//'cov1_m.diag'
    CALL bmcm(kanal,lgauss) ! (A^TC_d^-1A + C_m^-1)^-1   -> cov_m_dc
    if (errnr.ne.0) RETURN
    

    IF (lres) THEN

       fetxt = ramd(1:lnramd)//slash(1:1)//'res_m.diag'

       WRITE(*,'(a)')ACHAR(13)//'calculating RES = '//&
            '(A^TC_d^-1A + C_m^-1)^-1A^TC_d^-1A'
       WRITE(fprun,'(a)')'RES = (A^TC_d^-1A + C_m^-1)^-1A^TC_d^-1A'

       CALL bres(kanal)
       if (errnr.ne.0) RETURN

       IF (lcov2) THEN

          fetxt = ramd(1:lnramd)//slash(1:1)//'cov2_m.diag'

          WRITE(*,'(a)')ACHAR(13)//'calculating MCM_2 = '//&
               '(A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A (A^TC_d^-1A + C_m^-1)^-1'
          WRITE(fprun,'(a)')'MCM_2 = (A^TC_d^-1A + C_m^-1)'//&
               '^-1 A^TC_d^-1A (A^TC_d^-1A + C_m^-1)^-1'

          CALL bmcm2(kanal)
          if (errnr.ne.0) RETURN

       END IF              ! lcov2
    END IF                 ! lres

!!! $ may be we ant to write out not only main diagonals before freeing

    IF (ALLOCATED (ata)) DEALLOCATE (ata)
    IF (ALLOCATED (ata_reg)) DEALLOCATE (ata_reg)
    IF (ALLOCATED (cov_m)) DEALLOCATE (cov_m)
    
    fetxt = 'uncertainty calculation time'
    CALL TOC(c1,fetxt)

  END SUBROUTINE buncert

  SUBROUTINE bata(kanal)
!!!$    
!!!$    Unterprogramm berechnet A^HC_d^-1A
!!!$
!!!$    Copyright Andreas Kemna 2009
!!!$    Created by  Roland Martin                            02-Nov-2009
!!!$    
!!!$    Last changed      RM                                   20-Feb-2011
!!!$
!!!$.....................................................................
!!!$   PROGRAMMINTERNE PARAMETER:
!!!$   Hilfsvariablen 
    INTEGER                       :: kanal
    INTEGER                       :: i,j,k
    REAL,DIMENSION(:),ALLOCATABLE :: dig
    REAL                          :: dig_min,dig_max
    INTEGER                                      :: c1
    CHARACTER(80)                                :: csz
!!!$....................................................................

!!!$  A^TC_d^-1A

    errnr = 1
    OPEN (kanal,file=TRIM(fetxt),status='replace',err=999)
    errnr = 4

    ALLOCATE(dig(manz))

    ata = 0d0

    CALL TIC(c1)

    IF (ldc) THEN

       !$OMP PARALLEL DEFAULT (none) &
       !$OMP SHARED (ata,sensdc,wmatd,wdfak,dig,manz,nanz) &
       !$OMP PRIVATE (i,j,k)
       !$OMP DO SCHEDULE (GUIDED,CHUNK_0)
       DO k=1,manz
          !       write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//&
          !            'ATC_d^-1A/ ',REAL( k * (100./manz)),'%'
          DO j=k,manz ! fills upper triangle (k,j)
             DO i=1,nanz
                ata(k,j) = ata(k,j) + sensdc(i,k) * & 
                     sensdc(i,j) * wmatd(i) * DBLE(wdfak(i))
             END DO
             ata(j,k) = ata(k,j) ! fills lower triangle (k,j)
          END DO
          dig(k) = ata(k,k)
       END DO
       !$OMP END PARALLEL

    ELSE
       print*,'NON DC'
       !$OMP PARALLEL DEFAULT (none) &
       !$OMP SHARED (ata,sens,wmatd,wdfak,dig,manz,nanz) &
       !$OMP PRIVATE (i,j,k)
       !$OMP DO SCHEDULE (GUIDED,CHUNK_0)
       DO k=1,manz
          !       write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//&
          !            'ATC_d^-1A/ ',REAL( k * (100./manz)),'%'
          DO j=k,manz ! fills upper triangle (k,j)
             DO i=1,nanz
                ata(k,j) = ata(k,j) + DCONJG(sens(i,k)) * &
                     sens(i,j) * wmatd(i) * DBLE(wdfak(i))
             END DO
             ata(j,k) = ata(k,j) ! fills lower triangle (k,j)
          END DO
          dig(k) = ata(k,k)
       END DO
       !$OMP END PARALLEL
    END IF

    csz = 'ATA time'
    CALL TOC(c1,csz)

!!!$    get min/max
    dig_min = MINVAL(dig)
    dig_max = MAXVAL(dig)
!!!$    write out 
    WRITE (kanal,*)manz
    DO i=1,manz
       WRITE (kanal,*)dig(i),LOG10(dig(i))-LOG10(dig_max)
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
!!!$    Unterprogramm berechnet A^HC_d^-1A+lam*C_m
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
    INTEGER                                      :: c1
    CHARACTER(80)                                :: csz
!!!$.....................................................................

!!!$  A^TC_d^-1A+lamC_m

    errnr = 1
    open(kanal,file=TRIM(fetxt),status='replace',err=999)
    errnr = 4

    ALLOCATE(dig(manz))

    CALL TIC (c1)

    IF (ltri == 0) THEN

       !$OMP PARALLEL DEFAULT (none) &
       !$OMP SHARED (dig,ata,ata_reg,smatm,manz,lam,nx) &
       !$OMP PRIVATE (i,j)
       !$OMP DO SCHEDULE (GUIDED,CHUNK_0)

       DO j=1,manz
          write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)// &
               'ATC_d^-1A+reg/ ',REAL( j * (100./manz)),'%'
          DO i=1,j            ! lower triangle

             IF (i == j) THEN
                ata_reg(i,j) = ata(i,j) + lam * smatm(i,1)

                IF (i+1 < manz) ata_reg(i+1,j) = ata(i+1,j) + &
                     lam * smatm(i+1,2) ! nebendiagonale in x richtung
                IF (i+nx < manz) ata_reg(i+nx,j) = ata(i+nx,j) + &
                     lam * smatm(i+nx,3) ! nebendiagonale in z richtung
             ELSE
                ata_reg(i,j) = ata(i,j) ! only aTa
             END IF

             ata_reg(j,i) = ata_reg(i,j) ! upper triangle 

          END DO

          dig(j) = ata_reg(j,j)
       END DO
       
       !$OMP END PARALLEL

    ELSE IF (ltri == 3.OR.ltri == 4) THEN

       ata_reg= ata

       DO i=1,manz
          dig (i) = ata(i,i) + lam * smatm(i,1)
          ata_reg(i,i) = dig(i)
       END DO

    ELSE IF (ltri == 1.OR.ltri == 2.OR.(ltri > 4 .AND. ltri < 15)) THEN

       !$OMP PARALLEL DEFAULT (none) &
       !$OMP SHARED (dig,ata,ata_reg,smatm,manz,lam,nachbar,smaxs) &
       !$OMP PRIVATE (i,j)
       !$OMP DO SCHEDULE (GUIDED,CHUNK_0)

       DO j=1,manz

          write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)// &
               'ATC_d^-1A+reg/ ',REAL( j * (100./manz)),'%'

          DO i=1,j

             IF (i==j) THEN   ! sparse C_m

                ata_reg(i,j) = ata(i,j) + lam * smatm(i,smaxs+1)

                DO k=1,nachbar(i,smaxs+1)

                   IF (nachbar(i,k) /= 0) ata_reg(nachbar(i,k),j) = &
                        ata(nachbar(i,k),j) + lam * smatm(i,k)

                END DO

             ELSE

                ata_reg(i,j) = ata(i,j)

             END IF

             ata_reg(j,i) = ata_reg(i,j) ! upper triangle 

          END DO

          dig(j) = ata_reg(j,j)

       END DO
       
       !$OMP END PARALLEL

    ELSE IF (ltri == 15) THEN

       ata_reg= ata + lam * smatm ! for full C_m..w

       DO j=1,manz
          dig(j) = ata_reg(j,j)
       END DO

    END IF

    csz = 'ATA+RTR time'
    CALL TOC(c1,csz)

    dig_min = MINVAL(dig)
    dig_max = MAXVAL(dig)

    WRITE (kanal,*)manz,lam
    DO i=1,manz
       WRITE (kanal,*)dig(i),LOG10(dig(i))
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

    ALLOCATE (dig(manz)) !dig is a working array and contains main diagonal


!!!$    get time
    CALL TIC(c1)

    cov_m = ata_reg

    IF (.NOT.PRESENT(ols).OR..NOT.ols) THEN !default

       WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//'Factorization...'

       CALL CHOLD(cov_m,dig,manz,errnr,lverb)
       IF (errnr /= 0) THEN
          fetxt='CHOLD mcm :: matrix not pos definite..'
          PRINT*,'Zeile(',abs(errnr),')'
          errnr = 108
          RETURN
       END IF

       csz = 'Factorization time'
       CALL TOC(c1,csz)


       CALL TIC(c1)
       WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//'Inverting...'

       CALL LINVD(cov_m,dig,manz,lverb)

       WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//'Filling lower Cov...'

       !$OMP PARALLEL DEFAULT (none) &
       !$OMP SHARED (cov_m,manz,lverb) &
       !$OMP PRIVATE (i,j)
       !$OMP DO SCHEDULE (GUIDED,CHUNK_0)

       DO i= 1 , manz

          IF (lverb) WRITE (*,'(A,t45,F6.2,A)',ADVANCE='no')&
               ACHAR(13)//'Filling lower/',REAL( i * (100./manz)),'%'

          DO j = 1 , i - 1

             cov_m(i,j) = cov_m(j,i)

          END DO

       END DO

       !$OMP END PARALLEL

       csz = 'Inversion time'
       CALL TOC(c1,csz)


    ELSE
       WRITE (*,'(a)',ADVANCE='no')'Inverting Matrix (Gauss elemination)'

       CALL gauss_dble(cov_m,manz,errnr)

       IF (errnr /= 0) THEN
          fetxt = 'error matrix inverse not found'
          PRINT*,'Zeile::',cov_m(abs(errnr),:)
          PRINT*,'Spalte::',cov_m(:,abs(errnr))
          errnr = 108
          RETURN
       END IF

       csz = 'Inversion time'
       CALL TOC(c1,csz)

    END IF

    IF (lverb) THEN

       ALLOCATE (work(manz,manz),STAT=errnr)
       IF (errnr/=0) THEN
          fetxt = 'Allocation problem WORK in bmcm'
          errnr = 97
          RETURN
       END IF

!!$    !$OMP PARALLEL DEFAULT (none) SHARED (cov_m,ata_reg,work)
       work = MATMUL(cov_m,ata_reg)
!!$    !$OMP END PARALLEL
       DO i=1,manz
          IF (ABS(work(i,i) - 1d0) > 0.1) PRINT*,'bad approximation at parameter'&
               ,i,work(i,i)
       END DO
       DEALLOCATE (work)

    END IF

    DO i=1,manz
       dig(i) = cov_m(i,i)
    END DO

    dig_min = MINVAL(dig)
    dig_max = MAXVAL(dig)

    errnr = 1
    OPEN (kanal,file=TRIM(fetxt),status='replace',err=999)
    errnr = 4

    WRITE (kanal,*)manz,lam
    DO i=1,manz
       WRITE (kanal,*)dig(i),SQRT(dig(i))*1d2
    END DO

    WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
    WRITE (*,*)'Max/Min:',dig_max,'/',dig_min

    CLOSE(kanal)

    IF (lverb_dat) THEN

       errnr = 1
       open(kanal,file=TRIM(fetxt)//'_full',status='replace',err=999)
       errnr = 4
       DO i = 1, manz
          WRITE (kanal,*)espx(i),espy(i),(cov_m(i,j),j=i,manz)
       END DO
       
       CLOSE (kanal)
    END IF

    DEALLOCATE (dig)

    errnr = 0
999 RETURN

  END SUBROUTINE bmcm

  SUBROUTINE bres(kanal)
!!!$    
!!!$    Unterprogramm berechnet Aufloesungsmatrix
!!!$    Fuer beliebige Triangulierung
!!!$    RES = (A^HC_d^-1A + C_m^-1)^-1 A^HC_d^-1A
!!!$    wobei (A^HC_d^-1A + C_m^-1) bereits invertiert wurde (cov_m)
!!!$
!!!$    Copyright Andreas Kemna 2009
!!!$
!!!$    Created by Roland Martin                            02-Nov-2009
!!!$    
!!!$    Last changed      RM                             20-Feb-2010
!!!$.....................................................................
!!!$   PROGRAMMINTERNE PARAMETER:
!!!$   Hilfsvariablen 
    INTEGER                                      :: i,j,kanal
    REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE     :: dig
    REAL(KIND(0D0))                              :: dig_min,dig_max,dum
    INTEGER                                      :: c1
    CHARACTER(80)                                :: csz
!!!$.....................................................................

!!!$  cal!!!$RES = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A

    errnr = 1
    open(kanal,file=TRIM(fetxt),status='replace',err=999)
    errnr = 4
!!!$    get time
    CALL TIC(c1)

    ata_reg = MATMUL(cov_m,ata) ! that's it...

    csz = 'MATMUL time'
    CALL TOC(c1,csz)

    ALLOCATE (dig(manz)) !prepare to write out main diagonal

    DO i=1,manz
       dig(i) = ata_reg(i,i)
    END DO

    dig_min = MINVAL(dig)
    dig_max = MAXVAL(dig)

    WRITE (kanal,*)manz,lam
    DO i=1,manz
       WRITE (kanal,*)dig(i),LOG10(dig(i))-LOG10(dig_max)
    END DO

    WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
    WRITE (*,*)'Max/Min:',dig_max,'/',dig_min

    CLOSE(kanal)

    IF (lverb_dat) THEN

       errnr = 1
       open(kanal,file=TRIM(fetxt)//'_full',status='replace',err=999)
       errnr = 4
       DO i = 1, manz
          WRITE (kanal,*)espx(i),espy(i),(ata_reg(i,j),j=i,manz)
       END DO
       
       CLOSE (kanal)
    END IF

    DEALLOCATE (dig)

    errnr = 0
999 RETURN

  END SUBROUTINE bres

  SUBROUTINE bmcm2(kanal)
!!!$    
!!!$    Unterprogramm berechnet Modellkovarianz nach Fehlerfortpflanzung
!!!$    MCM = (A^HC_d^-1A + C_m^-1)^-1 A^HC_d^-1A (A^HC_d^-1A + C_m^-1)^-1
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
    INTEGER                                      :: i
    REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE     :: dig
    REAL(KIND(0D0))                              :: dig_min,dig_max
    INTEGER                                      :: c1
    CHARACTER(80)                                :: csz
!!!$.....................................................................
!!!$   vorher wurde schon 
!!!$   ata_reg_dc= (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A berechnet
!!!$   cov_m_dc= (A^TC_d^-1A + C_m^-1)^-1


    errnr = 1
    open(kanal,file=TRIM(fetxt),status='replace',err=999)
    errnr = 4

!!!$    get time
    CALL TIC(c1)

    ata = MATMUL (ata_reg,cov_m)
    
    csz = 'MATMUL time'
    CALL TOC(c1,csz)

    ALLOCATE (dig(manz))

    DO i=1,manz
       dig(i) = ata(i,i)
    END DO

    dig_min = MINVAL(dig) 
    dig_max = MAXVAL(dig)

    WRITE (kanal,*)manz,lam
    DO i=1,manz
       WRITE (kanal,*)dig(i),SQRT(dig(i))*1d2
    END DO

    WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
    WRITE (*,*)'Max/Min:',dig_max,'/',dig_min

    CLOSE(kanal)

    DEALLOCATE (dig)

    errnr = 0

999 RETURN

  END SUBROUTINE bmcm2
END MODULE bmcm_mod

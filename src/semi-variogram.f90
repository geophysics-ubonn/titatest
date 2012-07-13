PROGRAM semi_variogram

  IMPLICIT none

  !!$c     smallest variogram distance and tolerance
  REAL(KIND(0D0)) :: lag_unit,lag_tol
  !!$c     sill value
  REAL(KIND(0D0)) :: sill
!!$c     number of equidistant lags (nlag) = INT(grid_max / grid_min)
  INTEGER :: nrows,rown1,rown2
!!$c Row number of the data from data file
  INTEGER :: my_trafo
!!$c should we transform data or not? 0 = lin, 1 = ln, 2 = log10
  INTEGER :: nlag,ic_nlag
!!$c     lag vector (nlag)
  REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: lag
!!$c     number of headtail pairs for each semivariogram N(lag)(nlag)
  INTEGER,DIMENSION(:),ALLOCATABLE :: ngam
!!$c     experimental semivariogram
!!$c     gam(lag)=1/N(lag)/2 * sum_k^N(lag) (tail - head)**2.
  REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: gam
!!$c     parameter mean and variance
  REAL(KIND(0D0)) :: par_mean,par_vari
!!$c     head and tail
  REAL(KIND(0D0)) :: head,tail
!!$c     h distance
  REAL(KIND(0D0)) :: h
!!$c     th = tail - head
  REAL(KIND(0D0)) :: th
!!!$ grid and parameer variable storage
  REAL(KIND(0d0)),DIMENSION(:),ALLOCATABLE :: grid,para
!!!$ input storage
  REAL(KIND(0d0)),DIMENSION(:,:),ALLOCATABLE :: inp_buff
  REAL (KIND(0D0)) :: PI


  CHARACTER(256) :: input_filename,data_filename,&
       output_filename,cbuff

  INTEGER :: i,j,k
  INTEGER :: nlines
  INTEGER :: ifp1

  LOGICAL :: oki

  PI = ACOS(-1d0)


1 FORMAT(3(G10.3,3X),I10)
2 FORMAT(a,I10,a,F10.3)

  CALL GET_UNIT(ifp1)

  input_filename = 'semi-variogram.inp'

  INQUIRE(FILE=TRIM(input_filename),EXIST=oki)

  IF (.NOT. oki) THEN
     PRINT*,'create default '//TRIM(input_filename)//&
          ' with input data'
     OPEN (ifp1,FILE=TRIM(input_filename),STATUS='replace')
     WRITE (ifp1,'(a)')'# filename of input data'
     WRITE (ifp1,'(a)')'field.dat'
     WRITE (ifp1,'(a)')'# numbe rof rows, row of variable 1, row of variable  2'
     WRITE (ifp1,'(a)')'2 1  2'
     WRITE (ifp1,'(a)')'# sill, number of lags'
     WRITE (ifp1,'(a)')'# if number of lags = 0, we use'
     WRITE (ifp1,'(a)')'# nlag = ANINT(MAXVAL(ABS(grid))/'//&
          'MINVAL(ABS(grid)))'
     WRITE (ifp1,'(a)')'0.3 0'
     WRITE (ifp1,'(a)')'# lag(i) = i * MINVAL(ABS(grid))*2d0, i=0,21'
     WRITE (ifp1,'(a)')'# data trafo: 0=lin, 1=ln, 2=log10'
     WRITE (ifp1,'(a)')'0'
     STOP 
  END IF

  my_trafo = 0

  OPEN (ifp1,FILE=TRIM(input_filename),STATUS='old')
  CALL READ_COMMENTS(ifp1)
  READ (ifp1,*,ERR=999,END=1000)data_filename
  CALL READ_COMMENTS(ifp1)
  READ (ifp1,*,ERR=999,END=1000)nrows,rown1,rown2
  CALL READ_COMMENTS(ifp1)
  READ (ifp1,*,ERR=999,END=1000)lag_unit,nlag
  CALL READ_COMMENTS(ifp1)
  READ (ifp1,*,ERR=999,END=1000)my_trafo

  CLOSE (ifp1)

  output_filename = TRIM(data_filename)//'.svario'

  INQUIRE(FILE=TRIM(data_filename),EXIST=oki)
  IF (.NOT. oki) THEN
     PRINT*,'Error missing field data '//TRIM(data_filename)
     STOP
  END IF

  OPEN (ifp1,FILE=TRIM(data_filename),STATUS='old')
!!!$ deduce number of data lines
!!!$ oki is true here, since it was set by inquire
  nlines = 0
  DO WHILE (oki)
     READ (ifp1,*,END=100,ERR=100)cbuff
     nlines = nlines + 1
     CYCLE
100  EXIT
  END DO

  IF (nlines == 0) THEN
     PRINT*,'there is no data in '//TRIM(data_filename)
     STOP
  ELSE
     PRINT*,'there are ',nlines,' data points'
  END IF

  ALLOCATE (grid(nlines),para(nlines),inp_buff(nlines,nrows))
  
  REWIND(ifp1)

  READ (ifp1,*)(inp_buff(i,:),i=1,nlines)
  
  CLOSE (ifp1)

  grid = inp_buff(:,rown1)
  para = inp_buff(:,rown2)


  SELECT CASE (my_trafo)

  CASE (2)
     PRINT*,'Log10 Trafo'
     para = LOG10(para)

  CASE (1)
     PRINT*,'Log Trafo'
     para = LOG(para)

  END SELECT

  grid = grid
!  PRINT*,'grid::',grid
!  PRINT*,'para::',para


  IF (nlag == 0) THEN
     IF (nlines > 1) THEN
        lag_unit = grid(1) - grid(2)
     ELSE
        lag_unit = grid(1)
     END IF
     nlag = ANINT(MAXVAL(ABS(grid)) / lag_unit) - 1
     PRINT*,'Using self definitions for nlags:',nlag
  END IF

  lag_tol = lag_unit / 2d0

!!$ get some memory and initialize
  ALLOCATE (lag(nlag),ngam(nlag),gam(nlag))
  ngam = 0; gam = 0d0

!!!$ define lag vector
  DO k = 1, nlag
     lag(k) = k * lag_unit
  END DO

  PRINT*,'LAG STAT::',REAL(lag_unit),nlag,REAL(lag_tol),&
       REAL(MAXVAL(lag))

  WRITE (*,'(/a/)',ADVANCE='no')'Calculating VARIOGRAM'

  par_mean = SUM(para) / nlines

  par_vari = 0d0

  DO i = 1, nlines
     
     tail = para(i)

     par_vari = par_vari + (tail - par_mean)**2d0

     IF (ABS(tail) > 1e20) CYCLE

     DO j = 1, nlines

        IF (i == j) CYCLE

        head = para(j)

        IF (ABS(head) > 1e20) CYCLE 

        h = ABS(grid(i) - grid(j))

        th = head - tail

        IF (th < EPSILON(th)) CYCLE

        DO k = 1, nlag

!!$           gam(k) = gam(k) + th * th
!!$
!!$           ngam(k) = ngam(k) + 1
           IF ( ABS( lag(k) - h ) < lag_tol ) THEN
              ngam(k) = ngam(k) + 1
!              gam(k) = gam(k) * DBLE(ngam(k) - 1) ! multiply with previous value which was the devisor
              gam(k) = gam(k) + DBLE(th * th) ! add correctly
!              gam(k) = gam(k) / DBLE(ngam(k)) ! devide all throug hactual collect -> mean value
!!$              PRINT*,REAL(h),k,ngam(k),REAL(gam(k)),REAL(th*th),&
!!$                   REAL(gam(k) * DBLE(ngam(k) ))
           END IF
           

        END DO

     END DO
  END DO

!! $normalize semi variogram
  ic_nlag = 0
  DO k = 1,nlag
     IF (ngam(k) > 0) THEN
        gam(k) = gam(k) / DBLE(ngam(k)) / 2d0
        ic_nlag = ic_nlag + 1
     END IF
  END DO

  par_vari = par_vari / (nlines -1) 
  par_vari = MAX( par_vari ,1d-5)


!!!$ OUTPUT
  
  OPEN (ifp1,FILE=TRIM(output_filename),STATUS='replace')
  WRITE (ifp1,2)'#   lag(h)'//ACHAR(9)//&
       'exp. semivariogram     model ## ',ic_nlag,&
       ' / parameter variance ',par_vari
  DO k=1,nlag
!     IF (ngam(k) > 0) WRITE (ifp1,1)lag(k),gam(k),ngam(k)
     WRITE (ifp1,*)REAL(lag(k)),REAL(gam(k)),ngam(k)
  END DO
  CLOSE (ifp1)

  OPEN (ifp1,FILE='semi-variogram.dat',STATUS='replace')
  WRITE (ifp1,'(2(F10.5,1x),I5)')(REAL(lag(k)),REAL(gam(k)),ngam(k),k=1,nlag)
  CLOSE (ifp1)

  OPEN (ifp1,FILE='semi-variogram.gnu',STATUS='replace')
  WRITE (ifp1,'(a)')'set st da l'
  WRITE (ifp1,'(a)')'set grid'
  WRITE (ifp1,'(a)')"set out 'semi_variogram.ps'"
  WRITE (ifp1,'(a)')'set term pos col enh 20'
  WRITE (ifp1,'(a)')'set pointsize 1.2'
  WRITE (ifp1,'(a)')'set key bot right Left samplen .3'
  WRITE (ifp1,'(a)')"set xlab offset 0,0.5 'Lag h/[m]'"
!!$  WRITE (ifp1,'(a,2(F10.2,a))')'set xrange [',MINVAL(ABS(grid)),':',&
!!$       MAXVAL(ABS(grid))/2d0,']'
  WRITE (ifp1,'(a)')"set ylab offset 2,0 'sv(h)=1/2N(h)"//&
       " {/Symbol S}_i(Z(m_i)-Z(m_i+h))^2, {/Symbol g}(h) /[SI]'"
  WRITE (ifp1,*)'p"'//TRIM(output_filename)//'"'// &
       ' u 1:2 w p lc 1 pt 7 ti "sv(h)",',par_vari,&
       ' w l lc 0 lw 4 lt 2 ti "Variance (va)"'
  CLOSE (ifp1)

  STOP

999 PRINT*,'error reading '//TRIM(input_filename)
  STOP 
1000 PRINT*,'End of input file'
  
END PROGRAM semi_variogram

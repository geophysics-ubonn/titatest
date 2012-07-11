PROGRAM semi_variogram

  IMPLICIT none

  !!$c     smallest variogram distance and tolerance
  REAL(KIND(0D0)) :: lag_unit,lag_tol
!!$c     number of equidistant lags (nlag) = INT(grid_max / grid_min)
  INTEGER :: nlag
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
!!!$
  REAL(KIND(0d0)),DIMENSION(:),ALLOCATABLE :: grid,para

  CHARACTER(256) :: filename,cbuff

  INTEGER :: i,j,k
  INTEGER :: nlines
  INTEGER :: ifp1

  LOGICAL :: oki

1 FORMAT(3(G10.3,3X),I10)
2 FORMAT(a,I10,a,F10.3)

  ifp1 = 11

  filename = 'semi-variogram.inp'

  INQUIRE(FILE=TRIM(filename),EXIST=oki)

  IF (.NOT. oki) THEN
     PRINT*,'please create '//TRIM(filename)//' with data'
     STOP 
  END IF

  

  OPEN (ifp1,FILE=TRIM(filename),STATUS='old')
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
     PRINT*,'there is no data in '//TRIM(filename)
     STOP
  ELSE
     PRINT*,'there are ',nlines,' data points'
  END IF

  ALLOCATE (grid(nlines),para(nlines))
  
  REWIND(ifp1)

  READ (ifp1,*)(grid(i),para(i),i=1,nlines)
  
  CLOSE (ifp1)

!  PRINT*,'grid::',grid
!  PRINT*,'para::',para
  
  lag_unit = MINVAL(ABS(grid))

  nlag = ANINT(MAXVAL(ABS(grid)) / lag_unit) / 2

  lag_tol = lag_unit / 2d0

!!$ get some memory and initialize
  ALLOCATE (lag(nlag),ngam(nlag),gam(nlag))
  ngam = 0; gam = 0d0



!!!$ define lag vector
  DO k = 1, nlag
     lag(k) = k * lag_unit
  END DO

  PRINT*,'LAG STAT::',lag_unit,nlag,lag_tol

  WRITE (*,'(/a/)',ADVANCE='no')'Calculating VARIOGRAM'

  par_mean = SUM(para) / nlines

  par_vari = 0d0

  DO i = 1, nlines
     
     tail = para(i)

     par_vari = par_vari + (tail - par_mean)**2d0

     DO j = 1, nlines

        IF (i == j) CYCLE

        head = para(j)

        h = ABS(grid(i) - grid(j))

        th = (tail - head)**2d0

        DO k = 1, nlag

           IF ( ABS( lag(k) - h ) < lag_tol ) THEN
              ngam(k) = ngam(k) + 1
              gam(k) = gam(k) + th
           END IF
           

        END DO

     END DO
  END DO

!! $normalize semi variogram
  DO k = 1,nlag
     IF (ngam(k) > 0) gam(k) = gam(k) /ngam(k) / 2d0
  END DO

  par_vari = par_vari / (nlines -1) 
  par_vari = MAX( par_vari ,1d-5)


!!!$ OUTPUT

  OPEN (ifp1,FILE='semi-variogram.dat',STATUS='replace')
  WRITE (ifp1,2)'#   lag(h)'//ACHAR(9)//&
       'exp. semivariogram     model ## ',nlag,&
       ' / parameter variance ',par_vari
  DO i=1,nlag
     WRITE (ifp1,1)lag(i),gam(i),ngam(i)
  END DO
  CLOSE (ifp1)

  OPEN (ifp1,FILE='variogram.gnu',STATUS='replace')
  WRITE (ifp1,'(a)')'set st da l'
  WRITE (ifp1,'(a)')'set grid'
  WRITE (ifp1,'(a)')"set out 'variograms.ps'"
  WRITE (ifp1,'(a)')'set term pos col enh 20'
  WRITE (ifp1,'(a)')'set pointsize 1.2'
  WRITE (ifp1,'(a)')'set key bot right Left samplen .3'
  WRITE (ifp1,'(a)')"set xlab offset 0,0.5 'Lag h/[m]'"
  WRITE (ifp1,'(a,2(F10.2,a))')'set xrange [',MINVAL(ABS(grid)),':',&
       MAXVAL(ABS(grid))/2d0,']'
  WRITE (ifp1,'(a)')"set ylab offset 2,0 'sv(h)=1/2N(h)"//&
       " {/Symbol S}_i(Z(m_i)-Z(m_i+h))^2, {/Symbol g}(h) /[SI]'"
  WRITE (ifp1,*)'p'// &
       '"semi-variogram.dat" u 1:2 w p lc 1 pt 7 ti "sv(h)",',&
       par_vari,' w l lc 0 lw 4 lt 2 ti "Variance (va)"'
  CLOSE (ifp1)


END PROGRAM semi_variogram

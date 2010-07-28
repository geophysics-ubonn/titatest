      SUBROUTINE bvariogram
c     
c     Unterprogramm zum Bestimmen von Variogrammen
c     
c     Copyright by Andreas Kemna
c     
c     Erste Version von Roland Martin                          19-Mar-2010
c     
c     Letzte Aenderung                                         23-Mar-2010
c     
c.....................................................................

      USE invmod         ! fuer par
      USE variomodel
      USE sigmamod       ! fuer sigma
      USE modelmod       ! fuer manz
      USE elemmod        ! fuer grid_min,grid_max,etc

      IMPLICIT none

      INCLUDE 'parmax.fin'      ! fuer die felddefinitionen in elem.fin
      INCLUDE 'konv.fin'        ! fuer alfx/alfz
      INCLUDE 'err.fin'

c     PROGRAMMINTERNE PARAMETER:-------------------------------------------
c     Indexvariablen
      INTEGER :: i,j,ik,jk,ifp
c     Maximale knotenanzahl nicht entarteter Elemente
      INTEGER :: smaxs
c     Schwerpunktskoordinaten der Flaechenelemente und gitterabstaende
      REAL(KIND(0D0)) :: spx1,spx2,spy1,spy2
c     th = Tail - Head; hx,hy,h distances in each direction
      REAL(KIND(0D0)) :: th,tail,head,hx,hy,h,mid_par
c     korrelation length for variogram models
      REAL(KIND(0D0)) :: Ix,Iy
c     smallest variogram distance and tolerance
      REAL(KIND(0D0)) :: lag_unit,lag_tol
c     smallest variogram distance and tolerance anisotrop
      REAL(KIND(0D0)) :: lag_unit_x,lag_tol_x
      REAL(KIND(0D0)) :: lag_unit_y,lag_tol_y
c     number of equidistant lags (nlag) = INT(grid_max / grid_min)
      INTEGER :: nlag,nlag_x,nlag_y
c     lag vector (nlag)
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: lag,lag_x,lag_y
c     number of headtail pairs for each semivariogram N(lag)(nlag)
      INTEGER,DIMENSION(:),ALLOCATABLE :: ngam_x,ngam_y,ngam
c     experimental semivariogram
c     gam(lag)=1/N(lag)/2 * sum_k^N(lag) (tail - head)**2.
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: gam_x,gam_y,gam
c     variogram model
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: mgam_x,mgam_y,mgam
      CHARACTER(80) :: tmgam,tmgam_x,tmgam_y,mti,tgam,tcov
! mti stores a string for variogram statistics, like korrelation length
! tgam stores the output string of get_vario
      CHARACTER (11) :: tg
c-----------------------------------------------------------------------

      errnr = 4

      smaxs = selanz(1)
c     define linear equidistant lag vector

      lag_unit = grid_min
      lag_unit_x = grid_minx
      lag_unit_y = grid_miny

      nlag = ANINT(grid_max / lag_unit) / 2
      nlag_x = ANINT(grid_maxx / lag_unit_x) / 2
      nlag_y = ANINT(grid_maxy / lag_unit_y) / 2

      lag_tol = lag_unit * .5
      lag_tol_x = lag_unit_x * .5
      lag_tol_y = lag_unit_y * .5

c     get memory
      ALLOCATE (lag(nlag),gam(nlag),ngam(nlag),
     1     lag_x(nlag_x),gam_x(nlag_x),ngam_x(nlag),
     1     lag_y(nlag_y),gam_y(nlag_y),ngam_y(nlag),
     1     mgam(nlag),mgam_x(nlag_x),mgam_y(nlag_y),STAT=errnr)
      
      IF (errnr/=0) THEN
         fetxt = 'Allocation problem in bvariogram'
         WRITE (*,'(/a/)')TRIM(fetxt)
         errnr = 97
         RETURN
      END IF

      ngam_x = 0;ngam_y = 0;ngam = 0
      gam_x = 0.;gam_y = 0.;gam = 0.

! gets the current variogram function parameters
      CALL get_vario(Ix,Iy,fetxt,0) 
! now prepare title string of gnuplot plot
      WRITE (mti,'(3(a,F3.1))')'Integral length a=',
     1     SQRT(Ix**2.+Iy**2.),', ax=',Ix,', ay=',Iy
c     for postscript 
      tg = '{/Symbol g}'
      WRITE (tgam,'(a)')tg//'(h)='//TRIM(fetxt) 
      IF (ltri == 15) THEN ! only meaningful for stochastical regu..
         CALL get_vario(Ix,Iy,fetxt,1) ! gets teh 
c     for postscript
         WRITE (tcov,'(a)')'C(h)='//TRIM(fetxt)
      END IF
c     compute synthetic variogram model
      par_vari = 1.
      DO i=1,nlag
         lag(i) = i*lag_unit
         h = lag(i)
         WRITE (tmgam,'(a)')tg//'(h)'
         mgam(i) = mvario(h,h,par_vari)
      END DO
c     anisotrop
      DO i=1,nlag_x
         lag_x(i) = i*lag_unit_x
         h = lag_x(i)
         WRITE (tmgam_x,'(a)')tg//'(hx)'
         mgam_x(i) = mvario(h,0D0,par_vari)
      END DO
      DO i=1,nlag_y
         lag_y(i) = i*lag_unit_y
         h = lag_y(i)
         WRITE (tmgam_y,'(a)')tg//'(hy)'
         mgam_y(i) = mvario(0D0,h,par_vari)
      END DO

      mid_par = SUM(LOG10(DBLE(sigma(1:elanz)))) / manz
      PRINT*,'sigma mean',mid_par
      par_vari = 0.

c     Experimentelles semi-variogram
      DO i=1,elanz

         WRITE (*,'(a,1X,F6.2,A)',ADVANCE='no')ACHAR(13)//'variogram/',
     1        REAL(i*(100./elanz)),'%'

         spx1=0.;spy1=0.
         DO ik=1,smaxs
            spx1 = spx1 + sx(snr(nrel(i,ik)))
            spy1 = spy1 + sy(snr(nrel(i,ik)))
         END DO
         spx1 = spx1/smaxs; spy1 = spy1/smaxs

         tail = LOG10(DBLE(sigma(i))) ! lin val
         
         par_vari = par_vari + (tail - mid_par)**2.

         DO j=1,elanz

            IF (i==j) CYCLE

            spx2=0.;spy2=0.
            DO jk=1,smaxs
               spx2 = spx2 + sx(snr(nrel(j,jk)))
               spy2 = spy2 + sy(snr(nrel(j,jk)))
            END DO
            spx2 = spx2/smaxs; spy2 = spy2/smaxs

            head = LOG10(DBLE(sigma(j)))

            hx = ABS(spx1 - spx2)
            hy = ABS(spy1 - spy2)
            h = SQRT(hx**2. + hy**2.)

            th = (tail - head)**2.

            DO ik = 1,nlag
c     lag - lag_tol < h < lag + lag_tol
               IF (ABS(lag(ik) - h) < lag_tol) THEN
                  ngam(ik) = ngam(ik) + 1
                  gam(ik) = gam(ik) + th
               END IF
            END DO
c     anisotrop lag vector y-dir
            IF (hx < lag_tol_x) THEN
               DO ik = 1,nlag_y
                  IF (ABS(lag_y(ik) - hy) < lag_tol_y) THEN
                     ngam_y(ik) = ngam_y(ik) + 1
                     gam_y(ik) = gam_y(ik) + th
                  END IF
               END DO
            END IF
c     anisotrop lag vec x-dir
            IF (hy < lag_tol_y) THEN
               DO ik = 1,nlag_x
                  IF (ABS(lag_x(ik) - hx) < lag_tol_x) THEN
                     ngam_x(ik) = ngam_x(ik) + 1
                     gam_x(ik) = gam_x(ik) + th
                  END IF
               END DO
            END IF

         END DO                 ! inner loop j=1,elanz
      END DO                    ! outer loop i=1,elanz

      DO i=1,nlag
         gam(i) = gam(i) / ngam(i) / 2.
      END DO
      DO i=1,nlag_x
         gam_x(i) = gam_x(i) / ngam_x(i) / 2.
      END DO
      DO i=1,nlag_y
         gam_y(i) = gam_y(i) / ngam_y(i) / 2.
      END DO

c     sets parameter variance..
      par_vari = MAX(par_vari / MAX(1,manz - 1), 1.e-5)

      mgam = par_vari * mgam
      mgam_x = par_vari * mgam_x
      mgam_y = par_vari * mgam_y

      CALL get_unit(ifp)
      OPEN (ifp,FILE='inv.variogram_x',STATUS='replace',ERR=999)
      WRITE (ifp,'(a,I10)')'#   lag(x-dir)'//ACHAR(9)//
     1     'anisotrop exp. semivariogram    model ##',nlag_x
      DO i=1,nlag_x
         WRITE (ifp,'(3(G10.3,3X),I10)',ERR=999)
     1        lag_x(i),gam_x(i),mgam_x(i),ngam_x(i)
      END DO
      CLOSE (ifp)
      OPEN (ifp,FILE='inv.variogram_y',STATUS='replace',ERR=999)
      WRITE (ifp,'(a,I10)')'#   lag(y-dir)'//ACHAR(9)//
     1     'anisotrop exp. semivariogram   model ##',nlag_y
      DO i=1,nlag_y
         WRITE (ifp,'(3(G10.3,3X),I10)',ERR=999)
     1        lag_y(i),gam_y(i),mgam_y(i),ngam_y(i)
      END DO
      CLOSE (ifp)
      OPEN (ifp,FILE='inv.variogram',STATUS='replace',ERR=999)
      WRITE (ifp,'(a,I10,a,G10.3)')'#   lag(h)'//ACHAR(9)//
     1     'exp. semivariogram     model ## ',nlag,
     1     ' / parameter variance ',10**par_vari
      DO i=1,nlag
         WRITE (ifp,'(3(G10.3,3X),I10)',ERR=999)
     1        lag(i),gam(i),mgam(i),ngam(i)
      END DO
      CLOSE (ifp)

      OPEN (ifp,FILE='variogram.gnu',STATUS='replace',ERR=999)
      WRITE (ifp,'(a)')'set st da l'
      WRITE (ifp,'(a)')'set grid'
      WRITE (ifp,'(a)')"set out 'variograms.eps'"
      WRITE (ifp,'(a)')'set term pos col sol enh 20'
      WRITE (ifp,'(a)')'set key bot right Left'
      IF (ltri == 15) THEN ! only meaningful for stochastical regu..
         WRITE (ifp,'(a)')'set tit "'//TRIM(mti)//'\n'//
     1        TRIM(tgam)//'\n'//TRIM(tcov)//'"'
      ELSE
         WRITE (ifp,'(a)')'set tit "'//TRIM(mti)//'\n'//
     1        TRIM(tgam)//'"'
      END IF
      WRITE (ifp,'(a)')"set xlab offset 0,0.5 'Lag h/[m]'"
      WRITE (ifp,'(a,2(F10.2,a))')
     1     'set xrange [',grid_min,':',grid_max/2.,']'
      WRITE (ifp,'(a)')"set ylab offset 2,0 'sv(h)=1/2N(h)"//
     1     " {/Symbol S}_i(Z(m_i+h)-Z(m_i))^2'"
      WRITE (ifp,'(a,F10.3,a)')"plot"//
     1     "'inv.variogram' u 1:2 w lp lc 1 lw 3 ti 'sv(h)',"//
     1     "'inv.variogram' u 1:3 w l lc 1 lw 3 ti '"//TRIM(tmgam)//
     1     "','inv.variogram_x' u 1:2 w lp lc 2 lw 3 ti 'sv(hx)',"//
     1     "'inv.variogram_x' u 1:3 w l lc 2 lw 3 ti '"//TRIM(tmgam_x)//
     1     "','inv.variogram_y' u 1:2 w lp lc 3 lw 3 ti 'sv(hy)',"//
     1     "'inv.variogram_y' u 1:3 w l lc 3 lw 3 ti '"//TRIM(tmgam_y)//
     1     "',",par_vari,
     1     " w l lc 0 lw 4 ti 'Variance (va)'"
      CLOSE (ifp)
      
      par_vari = 10**par_vari

      fetxt = ''
      CALL SYSTEM('which gnuplot > tmp.gnuplot')
      OPEN (ifp,FILE='tmp.gnuplot',STATUS='old',ERR=999)
      READ (ifp,'(a)',ERR=100)fetxt
      CLOSE (ifp)
      fetxt = TRIM(ADJUSTL(fetxt))//' < variogram.gnu'
      IF (fetxt /= '') CALL SYSTEM (TRIM(fetxt))
      
 100  DEALLOCATE (gam_x,gam_y,gam)
      DEALLOCATE (ngam_x,ngam_y,ngam)
      DEALLOCATE (lag,lag_x,lag_y)

      errnr = 0
 999  RETURN

      END 

      subroutine bsmatmsto
c     
c     Unterprogramm belegt die Kovarianzmatrix.   
c     Neue Regularisierungsmatrix (stoch. Kovarianzmatrix).
c     
c     Copyright by Andreas Kemna 2009
c     
c     Erste Version von A. August/R. Martin                    03-Apr-2009
c     
c     Letzte Aenderung   RM                                    23-Nov-2009
c     
c.....................................................................
      
      USE alloci
      USE tic_toc
      USE variomodel

      IMPLICIT none
      
      INCLUDE 'parmax.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'dat.fin'
      INCLUDE 'model.fin'
      INCLUDE 'konv.fin'
      INCLUDE 'inv.fin'
      INCLUDE 'err.fin'
c.....................................................................

c     !!!
c     Kovarianzmatrix
      REAL (KIND(0D0)), DIMENSION(:,:),ALLOCATABLE  :: CovTT 
      INTEGER,DIMENSION(:),ALLOCATABLE              :: IPIV
c     inverse Kovarianzmatrix -> smatm
      integer*4            :: ErrorFlag

c     PROGRAMMINTERNE PARAMETER:
c     Schwerpunktskoordinaten der Flaechenelemente ij
      REAL(KIND(0D0)) :: spx1,spx2,spy1,spy2,h,sd_el
c     Korrelation lengths and variance (var)
      REAL(KIND(0D0)) :: Ix,Iy,var
!     gibt es evtl schon eine inverse?
      logical              :: ex,exc         
c     Hilfsvariablen
      integer              :: i,j,l,smaxs,ifp,c1,c2,se,mi,st,ta
c     clearscreen
      CHARACTER(80)        :: csz

      CALL get_unit(ifp)
      CALL TIC()

      var = 1.

      DO i=1,79
         csz(i:i+1)=' '
      END DO

      CALL gvario (Ix,Iy,csz)   ! get korrelation length and vario model string
      
      WRITE (*,'(A80)')ACHAR(13)//TRIM(csz)

      WRITE (*,'(A,1X,F6.2,1X,A)')ACHAR(13)//
     1     'Speicher fuer model covariance: ',
     1     REAL ((manz**2*8.)/(1024.**3)),'GB'
      
      IF (.NOT.ALLOCATED (smatm))
     1     ALLOCATE (smatm(manz,manz),STAT=errnr)
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem smatm in bsmatmsto'
         errnr = 97
         RETURN
      END IF

      smaxs=MAXVAL(selanz)
      
c     Belege die Matrix
c     covTT=0

      smatm = var

      INQUIRE(FILE='tmp.smatmi',EXIST=ex) ! already an inverse c_m ?

      IF (ex) THEN

         WRITE (*,'(a)',ADVANCE='no')'checking tmp.smatmi'
         OPEN (ifp,FILE='tmp.smatmi',STATUS='old',
     1        ACCESS='sequential',FORM='unformatted')
         READ (ifp) i
         IF (i == manz) THEN
            PRINT*,'seems appropriate inverse smatm'
            READ (ifp) smatm
         END IF
         CLOSE (ifp)

      ELSE

         do i = 1 , manz
            WRITE (*,'(a,1X,F6.2,A)',ADVANCE='no')ACHAR(13)//'cov/',
     1           REAL(i*(100./manz)),'%'
            spx1 = 0.;spy1 = 0.
            do l=1,smaxs
               spx1 = spx1 + sx(snr(nrel(i,l)))
               spy1 = spy1 + sy(snr(nrel(i,l)))
            end do
            spx1 = spx1 / smaxs ! x- schwerpunkt
            spy1 = spy1 / smaxs ! y- schwerpunkt

            do j = i , manz     ! fills upper triangle

               spx2 = 0.; spy2 = 0.
               DO l=1,smaxs
                  spx2 = spx2 + sx(snr(nrel(j,l)))
                  spy2 = spy2 + sy(snr(nrel(j,l)))
               END DO
               spx2 = spx2 / smaxs ! x- schwerpunkt
               spy2 = spy2 / smaxs ! y- schwerpunkt
               
               h = SQRT(((spx1 - spx2) / Ix)**2. +
     1              ((spy1 - spy2) / Iy)**2.)

               smatm(i,j) = mcova(h,var)

               IF (i /= j) smatm(j,i) = smatm(i,j) ! lower triangle

            end do
         end do


         PRINT*,'bestimme nun C_m^-1'
c     Berechne nun die Inverse der Covarianzmatrix!!!

c$$$  PRINT*,'   Gauss elemination ... '
         CALL gauss_dble(smatm,manz,errorflag)
c$$$         IF (.NOT.ALLOCATED (covTT)) ALLOCATE (covTT(manz,manz))
c$$$         covtt = 0.
c$$$         PRINT*,'   Cholesky factorization (Schwarz)... '
c$$$         CALL chold(smatm,covTT,manz,errorflag)
         IF (errorflag/=0) THEN
            fetxt='there was something wrong..'
            PRINT*,'Zeile(',abs(errorflag),
     1           ')::',smatm(abs(errorflag),:)
            PRINT*,'Spalte::',smatm(:,abs(errorflag))
            errnr = 108
         END IF

c$$$         PRINT*,'   Invertiere smatm ... '
c$$$         CALL linv(covTT,smatm,manz)

         IF (errorflag==0) THEN
            WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//
     1           'got inverse and write out'
            OPEN (ifp,FILE='tmp.smatmi',STATUS='replace',
     1           ACCESS='sequential',FORM='unformatted')
            WRITE (ifp) manz
            WRITE (ifp) smatm
            CLOSE (ifp)
         ELSE
            PRINT*,'got NO inverse'
            errnr = 108
            RETURN
         END IF     
      END IF

      IF (ALLOCATED (covTT)) DEALLOCATE (covTT)
      
      CALL TOC()

      END

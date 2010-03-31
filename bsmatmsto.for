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
c     smatm file name
      CHARACTER(10)        :: fsmat
c     clearscreen
      CHARACTER(80)        :: csz

      CALL get_unit(ifp)
      CALL TIC()

      var = 1.

      fsmat = 'inv.smatmi'

      DO i=1,79
         csz(i:i+1)=' '
      END DO

      CALL get_vario (Ix,Iy,csz,1)   
c!!!  get korrelation length and vario model string
      
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

      smatm = var

      INQUIRE(FILE=fsmat,EXIST=ex) ! already an inverse c_m ?

      IF (ex) THEN

         WRITE (*,'(a)',ADVANCE='no')'checking '//fsmat
         OPEN (ifp,FILE=fsmat,STATUS='old',
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

            do j = i+1 , manz     ! fills upper triangle

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

               smatm(j,i) = smatm(i,j) ! upper triangle

            end do
         end do


         PRINT*,'bestimme nun C_m^-1'
c     Berechne nun die Inverse der Covarianzmatrix!!!
         IF (lgauss) THEN
            PRINT*,'   Gauss elemination ... '
            CALL gauss_dble(smatm,manz,errorflag)
            IF (errorflag/=0) THEN
               fetxt='there was something wrong..'
               PRINT*,'Zeile(',abs(errorflag),
     1              ')::',smatm(abs(errorflag),:)
               PRINT*,'Spalte::',smatm(:,abs(errorflag))
               errnr = 108
            END IF
         ELSE ! default..
            WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//
     1           'Factorization...'
            CALL DPOTRF('U',manz,smatm,manz,errorflag)
            IF (errorflag/=0) THEN
               fetxt='there was something wrong..'
               PRINT*,'Zeile(',abs(errorflag),
     1              ')::',smatm(abs(errorflag),:)
               PRINT*,'Spalte::',smatm(:,abs(errorflag))
               errnr = 108
            END IF
            WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//
     1           'Inverting...'
            CALL DPOTRI('U',manz,smatm,manz,errnr)
            IF (errorflag/=0) THEN
               fetxt='there was something wrong..'
               PRINT*,'Zeile(',abs(errorflag),
     1              ')::',smatm(abs(errorflag),:)
               PRINT*,'Spalte::',smatm(:,abs(errorflag))
               errnr = 108
            END IF
            WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//
     1           'Filling lower C_m...'
            DO i= 1,manz
               WRITE (*,'(A,1X,F6.2,A)',ADVANCE='no')
     1              ACHAR(13)//ACHAR(9)//ACHAR(9)//
     1              ACHAR(9)//'/ ',REAL( i * (100./manz)),'%'
               DO j = i+1,manz
                  smatm(j,i)=smatm(i,j)
               END DO
            END DO
         END IF

         IF (errorflag==0) THEN
            WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//
     1           'got inverse and write out'
            OPEN (ifp,FILE=fsmat,STATUS='replace',
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

      CALL TOC()

      END

      subroutine bsmatmsto()
c
c     Unterprogramm belegt die Kovarianzmatrix.   
c     Neue Regularisierungsmatrix (stoch. Kovarianzmatrix).
c
c     Andreas Kemna                                            29-Feb-1996
c     Letzte Aenderung   03-Apr-2009

c.....................................................................
      
      USE alloci
      
      IMPLICIT none
      
      INCLUDE 'parmax.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'dat.fin'
      INCLUDE 'model.fin'
      INCLUDE 'konv.fin'
      INCLUDE 'inv.fin'

c.....................................................................

c     !!!
c     Varianz !!!!
      real :: var = 1

c     integral scale

      real :: Ix,Iz

c     Kovarianzmatrix
      REAL (KIND(0D0)), DIMENSION(:,:),ALLOCATABLE  :: CovTT 
      INTEGER,DIMENSION(:),ALLOCATABLE              :: IPIV
c     inverse Kovarianzmatrix -> smatm
      integer*4            :: ErrorFlag

c     PROGRAMMINTERNE PARAMETER:
c     Schwerpunktvektoren
c     !!
      real                * 8    xmeani, zmeani,xmeanj,zmeanj,
     1     h,sd_el,gamma
      
! gibt es evtl schon eine inverse?
      logical              :: ex,exc         
c     Hilfsvariablen
      integer              :: i,j,l,smaxs,ifp,c1,c2,se,mi,st,ta
c clearscreen
      CHARACTER(80)        :: csz

      CALL SYSTEM_CLOCK (c1,i)

      gamma = 2.

      DO i=1,79
         csz(i:i+1)=' '
      END DO

      IF (nz==1) THEN           !spherical
         WRITE (csz,'(a)')'Spherical model(va*(1- h*(1.5-.5*h**2)))'
      ELSE IF (nz==2) THEN      ! Gaussian
         WRITE (csz,'(a)')'Gaussian model(va*EXP(-3*h**2))'
      ELSE IF (nz==3) THEN      ! power
         WRITE (csz,'(a)')'Power model(va*dump**gamma)'
      ELSE                      ! exponential (default)
         WRITE (csz,'(a)')'Exponential model(va*EXP(-3*h))'
      END IF

      WRITE (*,'(A80)')ACHAR(13)//csz

      WRITE (*,'(A,1X,F6.2,1X,A)')ACHAR(13)//
     1     'Speicher fuer model covariance: ',
     1     REAL ((manz**2*8.)/(1024.**3)),'GB'
      
      IF (.NOT.ALLOCATED (smatm)) ALLOCATE (smatm(manz,manz))

      IF (alfx==0.) THEN
         Ix=esp_mit
         Iz=esp_mit
         PRINT*,'Choosing mean ESP distance as scale length:',Ix
      ELSE IF (alfz==0.) THEN
         Ix=esp_med
         Iz=esp_med
         PRINT*,'Choosing median ESP distance as scale length:',Ix
      ELSE
         Ix=alfx
         Iz=alfz
      END IF


      smaxs=MAXVAL(selanz)
      
c     Belege die Matrix
c     covTT=0

      smatm=0.
      do i = 1,manz
         WRITE (*,'(a,1X,F6.2,A)',ADVANCE='no')ACHAR(13)//'cov/',
     1        REAL(i*(100./manz)),'%'
         xmeani=0.
         do l=1,smaxs
            xmeani = xmeani + sx(snr(nrel(i,l)))
         end do
         xmeani = xmeani/smaxs  ! x- schwerpunkt

         zmeani = 0
         do l=1,smaxs
            zmeani = zmeani + sy(snr(nrel(i,l)))
         end do
         zmeani = zmeani/smaxs  ! y- schwerpunkt

         do j = 1,manz
            
            xmeanj = 0.
            do l=1,smaxs
               xmeanj = xmeanj + sx(snr(nrel(j,l)))
            end do
            xmeanj = xmeanj/smaxs ! x- schwerpunkt

            zmeanj = 0.
            do l=1,smaxs
               zmeanj = zmeanj + sy(snr(nrel(j,l)))
            end do
            zmeanj = zmeanj/smaxs ! y- schwerpunkt
            
            h = ((xmeani-xmeanj)/Ix)**2
     1           + ((zmeani-zmeanj)/Iz)**2
            
            h=sqrt(h)
            
            IF (nz==1) THEN     !spherical
               smatm(i,j) = var*(1. - h*(1.5 - .5*h**2.))
            ELSE IF (nz==2) THEN ! Gaussian
               smatm(i,j) = var*EXP(-3.*h**2.)
            ELSE IF (nz==3) THEN ! power
               smatm(i,j) = var*h**gamma
            ELSE                ! exponential (default)
               smatm(i,j) = var*EXP(-3.*h)
            END IF
         end do
      end do

c     Berechne nun (komponentenweise)
c     CovTT = sqrt(CovTT)

c     Berechne nun (komponentenweise)
c     CovTT = var*exp(-CovTT)


      exc=.TRUE.
      INQUIRE(FILE='tmp.smatmi',EXIST=ex)

      IF (ex) THEN
         PRINT*,'found tmp.smatmi'
         CALL get_unit(ifp)
         OPEN (ifp,FILE='tmp.smatmi',STATUS='old',
     1        ACCESS='sequential',FORM='unformatted')
         READ (ifp) i
         IF (i==manz) THEN
            PRINT*,'found appropriate inverse smatm'
            READ (ifp) smatm
         END IF
         exc=.FALSE.
         CLOSE (ifp)
      ELSE
      END IF

      IF (exc) THEN
         PRINT*,'bestimme nun inv{C_m}'
         IF (nx==-1) THEN
            PRINT*,'   DGESV (LAPACK)... '
            IF (.NOT.ALLOCATED (covTT)) ALLOCATE (covTT(manz,manz))
            IF (.NOT.ALLOCATED (IPIV)) ALLOCATE (IPIV(manz))
            covTT=0.0
            DO i=1,manz
               covTT(i,i)=1.0
            END DO
c$$$            CALL DPOTRF('U',manz,smatm,manz,errorflag)
c$$$            IF (errorflag/=0) THEN
c$$$               PRINT*,'there was something wrong..',errorflag
c$$$               STOP
c$$$            END IF
            PRINT*,'   Invertiere smatm ... '
c$$$            CALL MDPOTRI('U',manz,smatm,manz,errorflag)
c$$$            
c$$$            CALL DPOTRS('U',manz,manz,smatm,manz,covTT,
c$$$     1           manz,errorflag)
            CALL DGESV(manz,manz,smatm,manz,IPIV,covTT,
     1           manz,errorflag)
            IF (errorflag/=0) THEN
               PRINT*,'there was something wrong..'
               PRINT*,'Zeile::',covTT(abs(errorflag),:)
               PRINT*,'Spalte::',covTT(:,abs(errorflag))
               STOP
            END IF
            smatm=covTT
         ELSE IF (nx==-2) THEN        
            IF (.NOT.ALLOCATED (covTT)) ALLOCATE (covTT(manz,manz))
            PRINT*,'   Cholesky factorization (Schwarz)... '
            CALL chold(smatm,covTT,manz,errorflag)
            IF (errorflag/=0) THEN
               PRINT*,'there was something wrong..',errorflag
               STOP
            END IF
            PRINT*,'   Invertiere smatm ... '
            CALL linv(covTT,smatm,manz)
         ELSE IF (nx==-3) THEN        
            PRINT*,'   Find inv ... '
            IF (.NOT.ALLOCATED (covTT)) ALLOCATE (covTT(manz,manz))
            covTT=smatm
            CALL findinv(CovTT,smatm,manz,ErrorFlag)
         ELSE
            PRINT*,'   Gauss elemination ... '
c     Berechne nun die Inverse der Covarianzmatrix!!!
            CALL gauss(manz,errorflag)
         END IF
         IF (errorflag==0) THEN
            WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//
     1           'got inverse and write out'
            CALL get_unit(ifp)
            OPEN (ifp,FILE='tmp.smatmi',STATUS='replace',
     1           ACCESS='sequential',FORM='unformatted')
            WRITE (ifp) manz
            WRITE (ifp) smatm
            CLOSE (ifp)
         ELSE
            PRINT*,'got NO inverse'
            STOP
         END IF     
      END IF

c$$$      PRINT*,'Erasing border cell influence..'
c$$$      DO i=1,manz
c$$$         IF (nachbar(i,0)/=smaxs) THEN
c$$$            PRINT*,'Damping for cell#',i
c$$$            smatm(i,:)=0.
c$$$            smatm(i,i)=1.
c$$$         END IF
c$$$      END DO

      IF (ALLOCATED (covTT)) DEALLOCATE (covTT)

      CALL SYSTEM_CLOCK (c2,i)

      l = (c2-c1)/i ! Gesamt Sekunden
      mi=INT(l/60) ! Minuten
      st=INT(l/60/60) ! Stunden
      ta=INT(l/60/60/24) ! Tage
      se=l-mi*60-st*60*60-ta*60*60*24 ! Sekunden

 110  FORMAT(I3,'d/',1X,I2,'h/',1X,I2,'m/',1X,I2,'s')
      WRITE (*,110)ta,st,mi,se

      end

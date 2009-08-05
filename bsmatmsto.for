      subroutine bsmatmsto()

c     Unterprogramm belegt die Kovarianzmatrix.   Neue Regularisierungsmatrix (stoch. Kovarianzmatrix).

c     Andreas Kemna                                            29-Feb-1996
c     Letzte Aenderung   03-Apr-2009

c.....................................................................
c     !!
      USE alloci

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
      real*8,dimension(:,:),ALLOCATABLE  :: CovTT 
      
c     inverse Kovarianzmatrix -> smatm
      integer*4 :: ErrorFlag

c     PROGRAMMINTERNE PARAMETER:
c     Schwerpunktvektoren
c     !!
      real                * 8    xmeani, zmeani,xmeanj,zmeanj,
     1     dump

      logical :: ex,exc         ! gibt es evtl schon eine inverse?

c     Hilfsvariablen
      integer :: i,j,l,smaxs,ifp,c1,c2

      CALL SYSTEM_CLOCK (c1,i)

      IF (.NOT.ALLOCATED (smatm)) ALLOCATE (smatm(manz,manz))

      Ix=alfx
      Iz=alfz 

      smaxs=MAXVAL(selanz)

c     Belege die Matrix
c     covTT=0

      smatm=0.

      do i = 1,manz
         xmeani=0.
         do l=1,smaxs
            xmeani = xmeani + sx(snr(nrel(i,l)))
         end do
         xmeani = xmeani/smaxs  ! x- schwerpunkt

         zmeani=0
         do l=1,smaxs
            zmeani = zmeani + sy(snr(nrel(i,l)))
         end do
         zmeani = zmeani/smaxs  ! y- schwerpunkt

         do j = 1,manz
            
            xmeanj=0
            do l=1,smaxs
               xmeanj = xmeanj + sx(snr(nrel(j,l)))
            end do
            xmeanj = xmeanj/smaxs ! x- schwerpunkt

            zmeanj=0
            do l=1,smaxs
               zmeanj = zmeanj + sy(snr(nrel(j,l)))
            end do
            zmeanj = zmeanj/smaxs ! y- schwerpunkt
            
            dump = ((xmeani-xmeanj)/Ix)**2
     1           + ((zmeani-zmeanj)/Iz)**2
            
            dump=sqrt(dump)
            
            smatm(i,j) = var*EXP(-dump) 

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
      END IF

      IF (exc) THEN
         IF (nx<0) THEN
            PRINT*,''
            PRINT*,'   Cholesky factorization ... '
            PRINT*,''
            CALL DPOTRF('U',manz,smatm,manz,errorflag)
            IF (errorflag/=0) THEN
               PRINT*,'there was something wrong..',errorflag
               PRINT*,'Zeile::',smatm(abs(errorflag),:)
               PRINT*,'Spalte::',smatm(:,abs(errorflag))
               STOP
            END IF
            PRINT*,''
            PRINT*,'   Invertiere smatm ... '
            PRINT*,''
            CALL DPOTRI('U',manz,smatm,manz,errorflag)
            IF (errorflag/=0) THEN
               PRINT*,'there was something wrong..'
               PRINT*,'Zeile::',smatm(abs(errorflag),:)
               PRINT*,'Spalte::',smatm(:,abs(errorflag))
               STOP
            END IF
         ELSE IF (nz<0) THEN        
            PRINT*,''
            PRINT*,'   Find inv ... '
            PRINT*,''
            IF (.NOT.ALLOCATED (covTT)) ALLOCATE (covTT(manz,manz))
            covTT=smatm
            CALL findinv(CovTT,smatm,manz,ErrorFlag)
         ELSE
            PRINT*,''
            PRINT*,'   Gauss elemination ... '
            PRINT*,''
c     Berechne nun die Inierse der Covarianzmatrix!!!
            CALL gauss(smatm,manz,errorflag)
         END IF
         IF (errorflag==0) THEN
            PRINT*,'got inverse and write out'
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

      CALL SYSTEM_CLOCK (c2,i)
      WRITE (*,'(a,I6,a)')' in ',((c2-c1)/(i)),' s'

      end

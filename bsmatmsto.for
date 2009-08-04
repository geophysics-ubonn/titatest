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
      real                * 8    xmeani, zmeani,xmeanj,zmeanj

      logical :: ex,exc ! gibt es evtl schon eine inverse?

c     Hilfsvariablen
      integer :: i,j,l,smaxs,ifp
      

      IF (.NOT.ALLOCATED (covTT)) ALLOCATE (covTT(manz,manz))
      IF (.NOT.ALLOCATED (smatm)) ALLOCATE (smatm(manz,manz))

      Ix=alfx
      Iz=alfz 

      smaxs=MAXVAL(selanz)

c     Belege die Matrix
      covTT=0

      do i = 1,manz
         WRITE (*,'(A,I7)',ADVANCE='no')
     1        ACHAR(13)//ACHAR(9)//'Element ',i

         do l=1,smaxs
            xmeani = xmeani + sx(snr(nrel(i,l)))
         end do

         xmeani = xmeani/smaxs ! x- schwerpunkt

         do l=1,smaxs
            zmeani = zmeani + sy(snr(nrel(i,l)))
         end do

         zmeani = zmeani/smaxs ! y- schwerpunkt

         do j = 1,manz
            
            do l=1,smaxs
               xmeanj = xmeanj + sx(snr(nrel(j,l)))
            end do

            xmeanj = xmeanj/smaxs ! x- schwerpunkt

            do l=1,smaxs
               zmeanj = zmeanj + sy(snr(nrel(j,l)))
            end do
            
            zmeanj = zmeanj/smaxs ! y- schwerpunkt

c$$$            IF ((abs((xmeani-xmeanj)/xmeani)<0.3).OR.
c$$$     1           abs((zmeani-zmeanj)/zmeani)<0.3) CYCLE

            CovTT(i,j) = ((xmeani-xmeanj)/Ix)**2
     1           + ((zmeani-zmeanj)/Iz)**2

         end do
      end do


c     Berechne njn (komponentenweise)
      CovTT = sqrt(CovTT)

c     Berechne njn (komponentenweise)
      CovTT = var*exp(-CovTT)
      exc=.FALSE.
      INQUIRE(FILE='tmp.smatmi',EXIST=ex)

      IF (ex) THEN
         CALL get_unit(ifp)
         OPEN (ifp,FILE='tmp.smatmi',STATUS='old',
     1        ACCESS='sequential',FORM='unformatted')
         READ (ifp) i
         IF (i==manz) THEN
            PRINT*,'found appropriate inverse smatm'
            READ (ifp) smatm
         END IF
         CLOSE (ifp)
      END IF

      IF (exc) THEN
         WRITE (*,'(A)',ADVANCE='no')'   Invertiere CovTT ... '
c     Berechne nun die Inierse der Covarianzmatrix!!!
         CALL findinv(CovTT,smatm,manz,ErrorFlag)
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

      IF (ALLOCATED (covTT)) DEALLOCATE (covTT)

      return
      end

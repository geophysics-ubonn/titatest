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

      
c     Hilfsvariablen
      integer :: i,j,l,smaxi
      

      IF (.NOT.ALLOCATED (covTT)) ALLOCATE (covTT(manz,manz))
      IF (.NOT.ALLOCATED (smatm)) ALLOCATE (smatm(manz,manz))

      Ix=alfx
      Iz=alfz 

      smaxi=MAXVAL(selanz)

c     Belege die Matrix

      do i = 1,manz
         WRITE (*,'(A,I7)',ADVANCE='no')
     1        ACHAR(13)//ACHAR(9)//'Element ',i

         do l=1,smaxi
            xmeani = xmeani + sx(snr(nrel(i,l)))
         end do

         xmeani = xmeani/smaxi ! x- schwerpunkt

         do l=1,smaxi
            zmeani = zmeani + sy(snr(nrel(i,l)))
         end do

         zmeani = zmeani/smaxi ! y- schwerpunkt

         do j = 1,manz

            do l=1,smaxi
               xmeanj = xmeanj + sx(snr(nrel(j,l)))
            end do

            xmeanj = xmeanj/smaxi ! x- schwerpunkt

            do l=1,smaxi
               zmeanj = zmeanj + sy(snr(nrel(j,l)))
            end do
            
            zmeanj = zmeanj/smaxi ! y- schwerpunkt

            CovTT(i,j) = ((xmeani-xmeanj)/Ix)**2
     1           + ((zmeani-zmeanj)/Iz)**2

         end do
      end do


c     Berechne njn (komponentenweise)
      CovTT = sqrt(CovTT)

c     Berechne njn (komponentenweise)
      CovTT = var*exp(-CovTT)

      WRITE (*,'(A)',ADVANCE='no')'Invertiere CovTT ... '

c     Berechne nun die Inierse der Covarianzmatrix!!!
      CALL findinv(CovTT,smatm,manz,ErrorFlag)
      
      IF (errorflag==0) THEN
         PRINT*,'got inverse'
      ELSE
         PRINT*,'got NO inverse'
      end if      

      IF (ALLOCATED (covTT)) DEALLOCATE (covTT)

      return
      end

      program CDGR
      include 'omp_lib.h'

c     Declare variable types
      integer :: i, j, x, m, n
      double precision, dimension(:,:), allocatable :: W
      double precision :: cc,ss,time
      integer np, me

c     Get the size of matrix to use
      write(*,*) 'What size of matrix do you wish to use?'
      write(*,*) 'Number of rows (m) ='
      read(*,*) m
      write(*,*) 'Number of columns (n) ='
      read(*,*) n
      write(*,*)'m= ',m,' and n = ',n,' thank you.'

      allocate(W(m,n))

      do i=1,m
         do j=1,n
            W(i,j)=1000*rand(i+j)
         end do
      end do

!      do i=1,m
!         write(*,*)(W(i,j),j=1,n)
!      end do

      time=dtime(timearray)

c     Show time step (i) that each element would be annihilated during 
c     to leave the matrix W upper triangular

!$OMP parallel private(x,cc,ss,me,i) shared(W,np,m,n)
      np = omp_get_num_threads()
      me = omp_get_thread_num()

      do i=1,m+n-2
c     Every node uses the same value of i but the j values
c     are shared out and can be preformed at the same time 

!$OMP do schedule(dynamic,1)
         do j=1,n
            x=m+2*j-i-1
c     make sure element W(x,j) is with-in the matrix W
            if (j .lt. x) then 
               if (x .le. m) then
                   call drotg(W(x-1,j), W(x,j), cc, ss)
                   W(x,j)=0d0
                   call drot(n-j,W(x-1,j+1:n),1,W(x,j+1:n),1,cc,ss)
!                   W(x,j)=i
               endif
            endif
         enddo
!$OMP end do
      enddo
!$OMP end parallel

      time=dtime(timearray)
      write(*,*)'CDGR with ',np
      write(*,*)'m=',m,'n=',n,'time=',time,

c  Print W if you want to see how it was annihilated
!      do i=1,m
!         write(*,*)(W(i,j),j=1,n)
!      end do

      deallocate(W)

      stop
      end



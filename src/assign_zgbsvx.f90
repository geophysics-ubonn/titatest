subroutine assign_zgbsvx(ab,i,j,ku,n,value)
! documentation: http://www.netlib.org/lapack/explore-html/dc/d50/zgbsvx_8f.html
! Assign matrix entries i,j to band storage for ZGBSVX
! On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)
complex(kind(0d0)),dimension(2*ku+1,n),intent(inout) :: ab

integer :: i,j,ku,kl,n,index_i
complex(kind(0d0)) :: value
kl = ku
if ((max0(1,j-ku).le.i).and.(i.le.min0(n,j+kl))) then
    index_i = ku+1+i-j
    ab(index_i,j) = ab(index_i,j) + value
end if
end subroutine

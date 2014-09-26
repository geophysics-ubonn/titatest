subroutine assign_zpbsvx(ab,i,j,kd,n,value)
! documentation: http://www.netlib.org/lapack/explore-html/db/d51/zpbsvx_8f.html
! Assign matrix entries i,j to band storage for ZGBSVX
! On entry, the upper or lower triangle of the Hermitian band
!          matrix A, stored in the first KD+1 rows of the array, except
!          if FACT = 'F' and EQUED = 'Y', then A must contain the
!          equilibrated matrix diag(S)*A*diag(S).  The j-th column of A
!          is stored in the j-th column of the array AB as follows:
!          if UPLO = 'U', AB(KD+1+i-j,j) = A(i,j) for max(1,j-KD)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(N,j+KD).
use alloci, only: prec
complex(kind(0d0)),dimension(kd+1,n),intent(inout) :: ab

integer :: i,j,ku,kl,n,index_i
complex(kind(0d0)) :: value
character :: uplo
uplo = 'U'
if ((max0(1,j-kd).le.i).and.(i.le.j)) then
    index_i = kd+1+i-j
    ab(index_i,j) = ab(index_i,j) + value
end if
end subroutine

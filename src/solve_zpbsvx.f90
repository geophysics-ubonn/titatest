subroutine solve_zpbsvx(ab,x,b)
use elemmod, only:sanz,mb
implicit none
! Solve the linear system A*x = b in band storage mode using the LAPACK
! routine ZGBSVX.

! documentation: http://www.netlib.org/lapack/explore-html/dc/d50/zgbsvx_8f.html
character                           :: FACT
character                           :: TRANS
integer                             :: N
integer                             :: KD
integer                             :: NRHS
complex*16, dimension(mb+1,sanz)  :: AB
integer                             :: LDAB
complex*16, dimension(mb+1,sanz)  :: AFB
integer                             :: LDAFB
integer, dimension(sanz)            :: IPIV
character                           :: EQUED
double precision, dimension(sanz)   :: R
double precision, dimension(sanz)   :: C
complex*16, dimension(sanz,1)       :: B
integer                             :: LDB
complex*16, dimension(sanz,1)       :: X
integer                             :: LDX
double precision                    :: RCOND
double precision, dimension(1)      :: FERR
double precision, dimension(1)      :: BERR
complex*16, dimension(2*sanz)       :: WORK
double precision, dimension(sanz)   :: RWORK
integer                             :: INFO 
character                           :: UPLO
double precision, dimension(sanz)   :: S

uplo = 'U'
fact    = 'n'
trans   = 'n'
n       = sanz
kd      = mb
nrhs    = 1
!ab from input (alrealy allocated)
ldab    = kd+1
ldafb   = kd+1
ldb = sanz
ldx = sanz
!do it
call zpbsvx(fact, uplo, n, kd, nrhs, ab, ldab,&
             afb, ldafb, equed, s, b, ldb, x, ldx, rcond,&
             ferr, berr, work, rwork, info)
if ((info.ne.0).and.(info.le.sanz)) then
    print*,'ZGBSVx solver info var:',info
    print*,'info <= N:  U(i,i) is exactly zero.  The factorization'
    print*,'                   has been completed, but the factor U is exactly'
    print*,'                   singular, so the solution and error bounds'
    print*,'                   could not be computed. RCOND = 0 is returned.'
end if
if (info.eq.(sanz+1)) then
    print*,'ZGBSVx solver info: solution singular to machine precision'
    print*,'Condition number',rcond
end if
end subroutine

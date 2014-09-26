subroutine solve_zgbsvx(ab,x,b)
use elemmod, only:sanz,mb
use alloci, only:prec,lverbose
implicit none
! Solve the linear system A*x = b in band storage mode using the LAPACK
! routine ZGBSVX.

! documentation: http://www.netlib.org/lapack/explore-html/dc/d50/zgbsvx_8f.html
character                           :: FACT
character                           :: TRANS
integer                             :: N
integer                             :: KL
integer                             :: KU
integer                             :: NRHS
complex(prec), dimension(2*mb+1,sanz)  :: AB
integer                             :: LDAB
complex(prec), dimension(3*mb+1,sanz)  :: AFB
integer                             :: LDAFB
integer, dimension(sanz)            :: IPIV
character                           :: EQUED
double precision, dimension(sanz)   :: R
double precision, dimension(sanz)   :: C
complex(prec), dimension(sanz,1)       :: B
integer                             :: LDB
complex(prec), dimension(sanz,1)       :: X
integer                             :: LDX
double precision                    :: RCOND
double precision, dimension(1)      :: FERR
double precision, dimension(1)      :: BERR
complex(prec), dimension(2*sanz)       :: WORK
double precision, dimension(sanz)   :: RWORK
integer                             :: INFO 

x = 0D0
fact    = 'e'
trans   = 'n'
n       = sanz
kl      = mb
ku      = mb
nrhs    = 1
!ab from input (alrealy allocated)
ldab    = kl+ku+1
ldafb   = 2*kl+ku+1
ldb = sanz
ldx = sanz
!do it
call zgbsvx(fact, trans, n, kl, ku, nrhs, ab, ldab,&
             afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond,&
             ferr, berr, work, rwork, info)
if (lverbose) print*,'zgbsvx solver infos',info,'condition',rcond
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

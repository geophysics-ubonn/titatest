MODULE cjgmod
!!!$ ------------------------------------------------------------------
!!$ Module for the data of Conjugate gradient to solve Linear systems
!!$ Copyright by Andreas Kemna   2010
!!$ Edited first by Roland Martin                          27-07-2010
!!$ -------------------------------------------------------------------
!!$ Right hand side (RHS) vector
 COMPLEX(KIND(0D0)),PUBLIC,DIMENSION(:),ALLOCATABLE  :: bvec
!!$ dc real valued version ov RHS
 REAL(KIND(0D0)),PUBLIC,DIMENSION(:),ALLOCATABLE     :: bvecdc
!!$ residual vector
 COMPLEX(KIND(0D0)),PUBLIC,DIMENSION(:),ALLOCATABLE  :: rvec
!!$ residual vector dc
 REAL(KIND(0D0)),PUBLIC,DIMENSION(:),ALLOCATABLE     :: rvecdc
!!$ intermediate vector (stores Ap)
 COMPLEX(KIND(0D0)),PUBLIC,DIMENSION(:),ALLOCATABLE  :: pvec
!!$ intermediate vector (stores Ap) dc
 REAL(KIND(0D0)),PUBLIC,DIMENSION(:),ALLOCATABLE     :: pvecdc
!!$ CG residuals
 REAL(KIND(0D0)),PUBLIC,DIMENSION(:),ALLOCATABLE     :: cgres
!!$ storage
 REAL(KIND(0D0)),PUBLIC,DIMENSION(:),ALLOCATABLE     :: cgres2
!!$ CG Epsilon
 REAL(KIND(0D0)),PUBLIC                              :: eps
!!$ maximum number of CG steps
 INTEGER (KIND = 4),PUBLIC                           :: ncgmax
!!$ actual number of CG steps..
 INTEGER (KIND = 4),PUBLIC                           :: ncg
!!$c assist vectors
 REAL(KIND(0D0)),ALLOCATABLE,DIMENSION(:)            :: apdc
 COMPLEX(KIND(0D0)),ALLOCATABLE,DIMENSION(:)         :: ap
END MODULE cjgmod

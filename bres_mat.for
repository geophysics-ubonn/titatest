      SUBROUTINE bres_mat
      
      USE alloci
      USE femmod
      USE datmod
      USE invmod
      USE modelmod
      USE elemmod
      
      IMPLICIT none
      
      INCLUDE 'konv.fin'


!     help variable  
      INTEGER :: i,j,k,l

      
      IF (.NOT.ALLOCATED(atadc)) ALLOCATE (atadc(manz,manz))

      
      atadc=MATMUL(TRANSPOSE(sens),sens)
      
      
      
      IF (.NOT.ALLOCATED(inv_atadc)) 
     1     ALLOCATE (inv_atadc(manz,manz))
      
      END SUBROUTINE bres_mat
